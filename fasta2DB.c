/*******************************************************************************************
 *
 *  Add .fasta files to a DB:
 *     Adds the given fasta files in the given order to <path>.db.  If the db does not exist
 *     then it is created.  All .fasta files added to a given data base must have the same
 *     header format and follow Pacbio's convention.  A file cannot be added twice and this
 *     is enforced.  The command either builds or appends to the .<path>.idx and .<path>.bps
 *     files, where the index file (.idx) contains information about each read and their offsets
 *     in the base-pair file (.bps) that holds the sequences where each base is compessed
 *     into 2-bits.  The two files are hidden by virtue of their names beginning with a '.'.
 *     <path>.db is effectively a stub file with given name that contains an ASCII listing
 *     of the files added to the DB and possibly the block partitioning for the DB if DBsplit
 *     has been called upon it.
 *
 *  Author:  Gene Myers
 *  Date  :  May 2013
 *  Modify:  DB upgrade: now *add to* or create a DB depending on whether it exists, read
 *             multiple .fasta files (no longer a stdin pipe).
 *  Date  :  April 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-v] <path:db> ( -f<file> | -i[<name>] | <input:fasta> ... )";

static char number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };

typedef struct
  { int    argc;
    char **argv;
    FILE  *input;
    int    count;
    char  *name;
  } File_Iterator;

File_Iterator *init_file_iterator(int argc, char **argv, FILE *input, int first)
{ File_Iterator *it;

  it = Malloc(sizeof(File_Iterator),"Allocating file iterator");
  if (it == NULL)
    return (NULL);
  it->argc  = argc;
  it->argv  = argv;
  it->input = input;
  if (input == NULL)
    it->count = first;
  else
    { it->count = 1;
      rewind(input);
    }
  return (it);
}

int next_file(File_Iterator *it)
{ static char nbuffer[MAX_NAME+8];

  if (it->input == NULL)
    { if (it->count >= it->argc)
        return (0);
      it->name = it->argv[it->count++];
    }
  else
    { char *eol;

      if (fgets(nbuffer,MAX_NAME+8,it->input) == NULL)
        { if (feof(it->input))
            return (0);
          fprintf(stderr,"%s: IO error reading line %d of -f file of names\n",Prog_Name,it->count);
          it->name = NULL;
          return (1);
        }
      if ((eol = index(nbuffer,'\n')) == NULL)
        { fprintf(stderr,"%s: Line %d in file list is longer than %d chars!\n",
                         Prog_Name,it->count,MAX_NAME+7);
          it->name = NULL;
          return (1);
        }
      *eol = '\0';
      it->count += 1;
      it->name  = nbuffer;
    }
  return (1);
}


int main(int argc, char *argv[])
{ FILE  *istub, *ostub;
  char  *dbname;
  char  *root, *pwd;

  FILE  *bases, *indx;
  int64  boff, ioff;

  int    ifiles, ofiles, ocells;
  char **flist;

  DAZZ_DB db;
  int     ureads;
  int64   offset;

  char   *PIPE;
  FILE   *IFILE;
  int     VERBOSE;

  //   Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("fasta2DB")

    IFILE = NULL;
    PIPE  = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'f':
            IFILE = fopen(argv[i]+2,"r");
            if (IFILE == NULL)
              { fprintf(stderr,"%s: Cannot open file of inputs '%s'\n",Prog_Name,argv[i]+2);
                exit (1);
              }
            break;
          case 'i':
            PIPE = argv[i]+2;
            if (PIPE[0] != '\0')
              { FILE *temp;

                temp = fopen(PIPE,"w");
                if (temp == NULL)
                  { fprintf(stderr,"%s: Cannot create -i name '%s'\n",Prog_Name,argv[i]+2);
                    exit (1);
                  }
                fclose(temp);
                if (unlink(PIPE) != 0)
                  fprintf(stderr,"%s: [WARNING] Could not delete temporary file %s\n",
                                 Prog_Name,PIPE);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (IFILE != NULL && PIPE != NULL)
      { fprintf(stderr,"%s: Cannot use both -f and -i together\n",Prog_Name);
        exit (1);
      }

    if ( (IFILE == NULL && PIPE == NULL && argc <= 2) || 
        ((IFILE != NULL || PIPE != NULL) && argc != 2))
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -f: import files listed 1/line in given file.\n");
        fprintf(stderr,"      -i: import data from stdin, use optiona name as data source.\n");
        fprintf(stderr,"        : otherwise, import sequence of specified files.\n");
        exit (1);
      }
  }

  //  Try to open DB file, if present then adding to DB, otherwise creating new DB.  Set up
  //  variables as follows:
  //    dbname = full name of db = <pwd>/<root>.db
  //    istub  = open db file (if adding) or NULL (if creating)
  //    ostub  = new image of db file (will overwrite old image at end)
  //    bases  = .bps file positioned for appending
  //    indx   = .idx file positioned for appending
  //    ureads = # of reads currently in db
  //    offset = offset in .bps at which to place next sequence
  //    ioff   = offset in .idx file to truncate to if command fails
  //    boff   = offset in .bps file to truncate to if command fails
  //    ifiles = # of .fasta files to add
  //    ofiles = # of .fasta files added so far
  //    ocells = # of SMRT cells already in db
  //    flist  = [0..ifiles+ocells] list of file names (root only) added to db so far

  { int i;

    root   = Root(argv[1],".db");
    pwd    = PathTo(argv[1]);
    dbname = Strdup(Catenate(pwd,"/",root,".db"),"Allocating db name");
    if (dbname == NULL)
      exit (1);

    if (PIPE != NULL)
      ifiles = 1;
    else if (IFILE == NULL)
      ifiles = argc-2;
    else
      { File_Iterator *ng;

        ifiles = 0;
        ng = init_file_iterator(argc,argv,IFILE,2);
        if (ng == NULL)
          exit (1);
        while (next_file(ng))
          { if (ng->name == NULL)
              exit (1);
            ifiles += 1;
          }
        free(ng);
      }

    bases = NULL;
    indx  = NULL;
    ostub = NULL;
    ioff  = 0;
    boff  = 0;

    istub = fopen(dbname,"r");
    if (istub == NULL)
      { ocells = 0;

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"w+");
        indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"w+");
        if (bases == NULL || indx == NULL)
          goto error;

        fwrite(&db,sizeof(DAZZ_DB),1,indx);

        ureads  = 0;
        offset  = 0;
      }
    else
      { if (fscanf(istub,DB_NFILE,&ocells) != 1)
          { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
            exit (1);
          }

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"r+");
        indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r+");
        if (bases == NULL || indx == NULL)
          exit (1);

        if (fread(&db,sizeof(DAZZ_DB),1,indx) != 1)
          { if (ferror(indx))
              fprintf(stderr,"%s: System error, read failed\n",Prog_Name);
            else
              fprintf(stderr,"%s: File %s.idx is corrupted\n",Prog_Name,root);
            exit (1);
          }
        if (fseeko(bases,0,SEEK_END) < 0)
          SYSTEM_READ_ERROR
        if (fseeko(indx, 0,SEEK_END) < 0)
          SYSTEM_READ_ERROR

        ureads = db.ureads;
        offset = ftello(bases);
        boff   = offset;
        ioff   = ftello(indx);
        if (boff < 0 || ioff < 0)
          SYSTEM_READ_ERROR
      }

    flist  = (char **) Malloc(sizeof(char *)*(ocells+ifiles),"Allocating file list");
    ostub  = Fopen(Catenate(pwd,"/",root,".dbx"),"w+");
    if (ostub == NULL || flist == NULL)
      goto error;

    if (fprintf(ostub,DB_NFILE,ocells+ifiles) < 0)   //  Will write again with correct value at end
      { fprintf(stderr,"%s: System error, write failed\n",Prog_Name);
        goto error;
      }
    ofiles = 0;
    for (i = 0; i < ocells; i++)
      { int  last;
        char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
          { if (ferror(istub))
              fprintf(stderr,"%s: System error, read failed\n",Prog_Name);
            else
              fprintf(stderr,"%s: File %s.db is corrupted\n",Prog_Name,root);
            goto error;
          }
        if (ofiles == 0 || strcmp(flist[ofiles-1],fname) != 0)
          if ((flist[ofiles++] = Strdup(fname,"Adding to file list")) == NULL)
            goto error;
        if (fprintf(ostub,DB_FDATA,last,fname,prolog) < 0)
          { fprintf(stderr,"%s: System error, write failed\n",Prog_Name);
            goto error;
          }
      }
  }

  { int            maxlen;
    int64          totlen, count[4];
    int64          rmax;
    int            pmax;
    DAZZ_READ     *prec;
    char          *read;
    int            c;
    File_Iterator *ng = NULL;

    //  Buffer for reads all in the same well

    pmax = 100;
    prec = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*pmax,"Allocating record buffer");
    if (prec == NULL)
      goto error;

    //  Buffer for accumulating .fasta sequence over multiple lines

    rmax  = MAX_NAME + 60000;
    read  = (char *) Malloc(rmax+1,"Allocating line buffer");
    if (read == NULL)
      goto error;

    totlen = 0;              //  total # of bases in new .fasta files
    maxlen = 0;              //  longest read in new .fasta files
    for (c = 0; c < 4; c++)  //  count of acgt in new .fasta files
      count[c] = 0;

    //  For each new input source do

    if (PIPE == NULL)
      { ng = init_file_iterator(argc,argv,IFILE,2);  //  Setup to read .fasta's
        if (ng == NULL)                              //    from command line or file
          goto error;
      }

    while (PIPE != NULL || next_file(ng))
      { FILE *input;
        char  prolog[MAX_NAME];
        char *path, *core;
        int   eof;

        //  Open it: <path>/<core>.fasta if file, stdin otherwise with core = PIPE or "stdout"

        if (PIPE == NULL)

          { if (ng->name == NULL) goto error;

            path  = PathTo(ng->name);
            core  = Root(ng->name,".fasta");
            if ((input = Fopen(Catenate(path,"/",core,".fasta"),"r")) == NULL)
              goto error;
            free(path);
          }

        else

          { if (PIPE[0] == '\0')
              core  = Strdup("stdout","Allocating file name");
            else
              core  = Strdup(PIPE,"Allocating file name");
            if (core == NULL)
              goto error;
            input = stdin;
          }

        //  Get the header of the first line.  If the file is empty skip.

        eof = (fgets(read,MAX_NAME,input) == NULL);
        if (eof || strlen(read) < 1)
          { fclose(input);
            if (PIPE != NULL)
              { fprintf(stderr,"Standard input is empty, terminating!\n");
                break;
              }
            fprintf(stderr,"Skipping '%s', file is empty!\n",core);
            free(core);
            continue;
          }

        //  Check that core is not too long and name is unique or last source if PIPE'd

        if (strlen(core) >= MAX_NAME)
          { fprintf(stderr,"%s: File name over %d chars: '%.200s'\n",
                           Prog_Name,MAX_NAME,core);
            goto error;
          }

        { int j;

 
          if (PIPE == NULL || (strcmp(core,"stdout") != 0 &&
                 (ofiles == 0 || strcmp(core,flist[ofiles-1]) != 0)))
            for (j = 0; j < ofiles; j++)
              if (strcmp(core,flist[j]) == 0)
                { fprintf(stderr,"%s: File %s.fasta is already in database %s.db\n",
                                 Prog_Name,core,Root(argv[1],".db"));
                  goto error;
                }
        }

        //   Add the file name to flist

        if (VERBOSE)
          { if (PIPE != NULL && PIPE[0] == '\0')
              fprintf(stderr,"Adding reads from stdio ...\n");
            else
              fprintf(stderr,"Adding '%s.fasta' ...\n",core);
            fflush(stderr);
          }
        flist[ofiles++] = core;

        // Check that the first line is a header and has PACBIO format.

        if (read[strlen(read)-1] != '\n')
          { fprintf(stderr,"File %s.fasta, Line 1: Fasta line is too long (> %d chars)\n",
                           core,MAX_NAME-2);
            goto error;
          }
        if (!eof && read[0] != '>')
          { fprintf(stderr,"File %s.fasta, Line 1: First header in fasta file is missing\n",core);
            goto error;
          }

        { char *find;
          int   well, beg, end, qv;

          find = index(read+1,'/');
          if (find != NULL && sscanf(find+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv) >= 3)
            { *find = '\0';
              strcpy(prolog,read+1);
              *find = '/';
            }
          else if (find != NULL && sscanf(find+1,"%d/ccs\n",&well) >= 1)
            { *find = '\0';
              strcpy(prolog,read+1);
              *find = '/';
            }
          else
            { fprintf(stderr,"File %s.fasta, Line 1: Pacbio header line format error\n",core);
              goto error;
            }
        }

        //  Read in all the sequences until end-of-file

        { int i, x;
          int nline, pwell, rlen, pcnt;

          pcnt  = 0;
          rlen  = 0;
          nline = 1;
          pwell = -1;
          while (!eof)
            { int   beg, end, clen;
              int   well, qv;
              char *find;

              find = index(read+(rlen+1),'/');
              if (find == NULL)
                { fprintf(stderr,"File %s.fasta, Line %d: Pacbio header line format error\n",
                                 core,nline);
                  goto error;
                }
              *find = '\0';
              if (strcmp(read+(rlen+1),prolog) != 0)
                { fprintf(ostub,DB_FDATA,ureads,core,prolog);
                  ocells += 1;
                  strcpy(prolog,read+(rlen+1));
                }
              *find = '/';
              x = sscanf(find+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv);
              if (x < 3)
                { char *secn = index(find+1,'/');
                  x = sscanf(find+1,"%d/ccs\n",&well);
                  if (secn == NULL || strncmp(secn+1,"ccs",3) != 0 || x < 1)
                    { fprintf(stderr,"File %s.fasta, Line %d: Pacbio header line format error\n",
                                     core,nline);
                      goto error;
                    }
                  beg = 0;
                  qv  = 0;
                }
              else if (x == 3)
                qv = 0;

              rlen  = 0;
              while (1)
                { eof = (fgets(read+rlen,MAX_NAME,input) == NULL);
                  nline += 1;
                  x = strlen(read+rlen)-1;
                  if (read[rlen+x] != '\n')
                    { if (read[rlen] == '>')
                        { fprintf(stderr,"File %s.fasta, Line %d:",core,nline);
                          fprintf(stderr," Fasta header line is too long (> %d chars)\n",
                                         MAX_NAME-2);
                          goto error;
                        }
                      else
                        x += 1;
                    }
                  if (eof || read[rlen] == '>')
                    break;
                  rlen += x;
                  if (rlen + MAX_NAME > rmax)
                    { rmax = ((int64) (1.2 * rmax)) + 1000 + MAX_NAME;
                      read = (char *) realloc(read,rmax+1);
                      if (read == NULL)
                        { fprintf(stderr,"File %s.fasta, Line %d:",core,nline);
                          fprintf(stderr," Out of memory (Allocating line buffer)\n");
                          goto error;
                        }
                    }
                }
              read[rlen] = '\0';

              for (i = 0; i < rlen; i++)
                { x = number[(int) read[i]];
                  count[x] += 1;
                  read[i]   = (char) x;
                }
              ureads += 1;
              totlen += rlen;
              if (rlen > maxlen)
                maxlen = rlen;

              prec[pcnt].origin = well;
              prec[pcnt].fpulse = beg;
              prec[pcnt].rlen   = rlen;
              prec[pcnt].boff   = offset;
              prec[pcnt].coff   = -1;
              prec[pcnt].flags  = qv;

              Compress_Read(rlen,read);
              clen = COMPRESSED_LEN(rlen);
              fwrite(read,1,clen,bases);
              offset += clen;

              if (pwell == well)
                { prec[pcnt].flags |= DB_CSS;
                  pcnt += 1;
                  if (pcnt >= pmax)
                    { pmax = ((int) (pcnt*1.2)) + 100;
                      prec = (DAZZ_READ *) realloc(prec,sizeof(DAZZ_READ)*pmax);
                      if (prec == NULL)
                        { fprintf(stderr,"File %s.fasta, Line %d: Out of memory",core,nline);
                          fprintf(stderr," (Allocating read records)\n");
                          goto error;
                        }
                    }
                }
              else if (pcnt == 0)
                pcnt += 1;
              else
                { x = 0;
                  for (i = 1; i < pcnt; i++)
                    if (prec[i].rlen > prec[x].rlen)
                      x = i;
                  prec[x].flags |= DB_BEST;
                  fwrite(prec,sizeof(DAZZ_READ),pcnt,indx);
                  prec[0] = prec[pcnt];
                  pcnt = 1;
                }
              pwell = well;
            }

          //  Complete processing of .fasta file: flush last well group, write file line
          //      in db image, and close file

          x = 0;
          for (i = 1; i < pcnt; i++)
            if (prec[i].rlen > prec[x].rlen)
              x = i;
          prec[x].flags |= DB_BEST;
          fwrite(prec,sizeof(DAZZ_READ),pcnt,indx);
        }

        fprintf(ostub,DB_FDATA,ureads,core,prolog);
        ocells += 1;

        if (input != stdin)
          fclose(input);
        else
          break;
      }

    //  Finished loading all sequences: update relevant fields in db record

    db.ureads = ureads;
    if (istub == NULL)
      { for (c = 0; c < 4; c++)
          db.freq[c] = (float) ((1.*count[c])/totlen);
        db.totlen = totlen;
        db.maxlen = maxlen;
        db.cutoff = -1;
        db.allarr = 0;
      }
    else
      { for (c = 0; c < 4; c++)
          db.freq[c] = (float) ((db.freq[c]*db.totlen + (1.*count[c]))/(db.totlen + totlen));
        db.totlen += totlen;
        if (maxlen > db.maxlen)
          db.maxlen = maxlen;
      }
  }

  //  If db has been previously partitioned then calculate additional partition points and
  //    write to new db file image

  if (db.cutoff >= 0)
    { int64      totlen, dbpos, size;
      int        nblock, ireads, tfirst, rlen;
      int        ufirst, cutoff, allflag;
      DAZZ_READ  record;
      int        i;

      if (VERBOSE)
        { fprintf(stderr,"Updating block partition ...\n");
          fflush(stderr);
        }

      //  Read the block portion of the existing db image getting the indices of the first
      //    read in the last block of the exisiting db as well as the partition parameters.
      //    Copy the old image block information to the new block information (except for
      //    the indices of the last partial block)

      if (fscanf(istub,DB_NBLOCK,&nblock) != 1)
        { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
          goto error;
        }
      dbpos = ftello(ostub);
      fprintf(ostub,DB_NBLOCK,0);
      if (fscanf(istub,DB_PARAMS,&size,&cutoff,&allflag) != 3)
        { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
          goto error;
        }
      fprintf(ostub,DB_PARAMS,size,cutoff,allflag);
      if (allflag)
        allflag = 0;
      else
        allflag = DB_BEST;

      nblock -= 1;
      for (i = 0; i <= nblock; i++)
        { if (fscanf(istub,DB_BDATA,&ufirst,&tfirst) != 2)
            { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
              goto error;
            }
          fprintf(ostub,DB_BDATA,ufirst,tfirst);
        }

      //  Seek the first record of the last block of the existing db in .idx, and then
      //    compute and record partition indices for the rest of the db from this point
      //    forward.

      fseeko(indx,sizeof(DAZZ_DB)+sizeof(DAZZ_READ)*ufirst,SEEK_SET);
      totlen = 0;
      ireads = 0;
      for (i = ufirst; i < ureads; i++)
        { if (fread(&record,sizeof(DAZZ_READ),1,indx) != 1)
            { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
              goto error;
            }
          rlen = record.rlen;
          if (rlen >= cutoff && (record.flags & DB_BEST) >= allflag)
            { ireads += 1;
              tfirst += 1;
              totlen += rlen;
              if (totlen >= size)
                { fprintf(ostub," %9d %9d\n",i+1,tfirst);
                  totlen = 0;
                  ireads = 0;
                  nblock += 1;
                }
            }
        }

      if (ireads > 0)
        { fprintf(ostub,DB_BDATA,ureads,tfirst);
          nblock += 1;
        }

      db.treads = tfirst;

      fseeko(ostub,dbpos,SEEK_SET);
      fprintf(ostub,DB_NBLOCK,nblock);    //  Rewind and record the new number of blocks
    }
  else
    db.treads = ureads;

  rewind(indx);
  fwrite(&db,sizeof(DAZZ_DB),1,indx);   //  Write the finalized db record into .idx

  rewind(ostub);                        //  Rewrite the number of files actually added
  fprintf(ostub,DB_NFILE,ocells);

  if (istub != NULL)
    fclose(istub);
  fclose(ostub);
  fclose(indx);
  fclose(bases);

  rename(Catenate(pwd,"/",root,".dbx"),dbname);   //  New image replaces old image

  exit (0);

  //  Error exit:  Either truncate or remove the .idx and .bps files as appropriate.
  //               Remove the new image file <pwd>/<root>.dbx

error:
  if (ioff != 0)
    { fseeko(indx,0,SEEK_SET);
      if (ftruncate(fileno(indx),ioff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.idx after error, truncate failed\n",
                       Prog_Name,root);
    }
  if (boff != 0)
    { fseeko(bases,0,SEEK_SET);
      if (ftruncate(fileno(bases),boff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.bps after error, truncate failed\n",
                       Prog_Name,root);
    }
  if (indx != NULL)
    { fclose(indx);
      if (ioff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".idx"));
    }
  if (bases != NULL)
    { fclose(bases);
      if (boff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".bps"));
    }
  if (ostub != NULL)
    { fclose(ostub);
      unlink(Catenate(pwd,"/",root,".dbx"));
    }
  if (istub != NULL)
    fclose(istub);

  exit (1);
}
