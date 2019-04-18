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

static char *Usage = "[-v] <path:dam> ( -f<file> | -i[<name>] | <input:fasta> ... )";

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
      0, 0, 0, 0, 0, 0, 4, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 4, 0,
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

  FILE  *bases, *indx, *hdrs;
  int64  boff, ioff, hoff, noff;

  int    ifiles, ofiles;
  char **flist;

  DAZZ_DB db;
  int     ureads;
  int64   offset, hdrset;

  char   *PIPE;
  FILE   *IFILE;
  int     VERBOSE;

  //   Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("fasta2DAM")

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
                unlink(PIPE);
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

  //  Try to open DAM file, if present then adding to DAM, otherwise creating new DAM.  Set up
  //  variables as follows:
  //    dbname = full name of map index = <pwd>/<root>.dam
  //    istub  = open db file (if adding) or NULL (if creating)
  //    ostub  = new image of db file (will overwrite old image at end)
  //    bases  = .bps file positioned for appending
  //    indx   = .idx file positioned for appending
  //    hdrs   = .hdr file positioned for appending
  //    ureads = # of reads currently in db
  //    offset = offset in .bps at which to place next sequence
  //    hdrset = offset in .hdr at which to place next header
  //    ioff   = offset in .idx file to truncate to if command fails
  //    boff   = offset in .bps file to truncate to if command fails
  //    hoff   = offset in .hdr file to truncate to if command fails
  //    ifiles = # of .fasta files to add
  //    ofiles = # of .fasta files added so far
  //    flist  = [0..ifiles+ofiles] list of file names (root only) added to dam so far

  { int i;

    root   = Root(argv[1],".dam");
    pwd    = PathTo(argv[1]);
    dbname = Strdup(Catenate(pwd,"/",root,".dam"),"Allocating map index name");
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
    hdrs  = NULL;
    ostub = NULL;
    ioff  = 0;
    boff  = 0;
    hoff  = 0;

    istub = fopen(dbname,"r");
    if (istub == NULL)
      { ofiles = 0;

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"w+");
        indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"w+");
        hdrs  = Fopen(Catenate(pwd,PATHSEP,root,".hdr"),"w+");
        if (bases == NULL || indx == NULL || hdrs == NULL)
          goto error;

        fwrite(&db,sizeof(DAZZ_DB),1,indx);

        ureads  = 0;
        offset  = 0;
        hdrset  = 0;
        boff    = 0;
        ioff    = 0;
        hoff    = 0;
      }
    else
      { if (fscanf(istub,DB_NFILE,&ofiles) != 1)
          { fprintf(stderr,"%s: %s.dam is corrupted, read failed 1\n",Prog_Name,root);
            goto error;
          }

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"r+");
        indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r+");
        hdrs  = Fopen(Catenate(pwd,PATHSEP,root,".hdr"),"r+");
        if (bases == NULL || indx == NULL || hdrs == NULL)
          goto error;

        if (fread(&db,sizeof(DAZZ_DB),1,indx) != 1)
          { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
            goto error;
          }
        fseeko(bases,0,SEEK_END);
        fseeko(indx, 0,SEEK_END);
        fseeko(hdrs, 0,SEEK_END);

        ureads = db.ureads;
        offset = ftello(bases);
        hdrset = ftello(hdrs);
        boff   = offset;
        ioff   = ftello(indx);
        hoff   = hdrset;
      }

    flist  = (char **) Malloc(sizeof(char *)*(ofiles+ifiles),"Allocating file list");
    ostub  = Fopen(Catenate(pwd,"/",root,".dbx"),"w+");
    if (ostub == NULL || flist == NULL)
      goto error;

    fprintf(ostub,DB_NFILE,ofiles+ifiles);
    noff = 0;
    for (i = 0; i < ofiles; i++)
      { int  last;
        char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
          { fprintf(stderr,"%s: %s.dam is corrupted, read failed 5(%d)\n",Prog_Name,root,i);
            goto error;
          }
        if ((flist[i] = Strdup(fname,"Adding to file list")) == NULL)
          goto error;
        noff = ftello(ostub);
        fprintf(ostub,DB_FDATA,last,fname,prolog);
      }
  }

  { int            maxlen;
    int64          totlen, count[4];
    int64          rmax;
    DAZZ_READ      prec;
    char          *read;
    int            append;
    int            c;
    File_Iterator *ng = NULL;

    //  Buffer for accumulating .fasta sequence over multiple lines

    rmax  = MAX_NAME + 60000;
    read  = (char *) Malloc(rmax+1,"Allocating line buffer");
    if (read == NULL)
      goto error;

    totlen = 0;              //  total # of bases in new .fasta files
    maxlen = 0;              //  longest read in new .fasta files
    for (c = 0; c < 4; c++)  //  count of acgt in new .fasta files
      count[c] = 0;

    //  For each .fasta file do:

    if (PIPE == NULL)
      { ng = init_file_iterator(argc,argv,IFILE,2);
        if (ng == NULL)
          goto error;
      }

    while (PIPE != NULL || next_file(ng))
      { FILE *input;
        char *path, *core;
        int   nline, eof, rlen;

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

        //  Check that core is not too long and name is unique or last source if PIPE'd
        //    If PIPE'd and last source, then overwrite last file line of new stub file.

        if (strlen(core) >= MAX_NAME)
          { fprintf(stderr,"%s: File name over %d chars: '%.200s'\n",
                           Prog_Name,MAX_NAME,core);
            goto error;
          }

        { int j;

          append = 0;
          if (PIPE == NULL || (strcmp(core,"stdout") != 0 &&
                 (ofiles == 0 || strcmp(core,flist[ofiles-1]) != 0)))
            { for (j = 0; j < ofiles; j++)
                if (strcmp(core,flist[j]) == 0)
                  { fprintf(stderr,"%s: File %s.fasta is already in database %s.dam\n",
                                   Prog_Name,core,Root(argv[1],".dam"));
                    goto error;
                  }
            }
	  else if (ofiles > 0 && strcmp(core,flist[ofiles-1]) == 0)
            { fseeko(ostub,noff,SEEK_SET);
              append = 1;
            }
        }

        //  Get the header of the first line.  If the file is empty skip.

        rlen  = 0;
        nline = 1;
        eof   = (fgets(read,MAX_NAME,input) == NULL);
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

        //   Add the file name to flist

        if (VERBOSE)
          { if (PIPE != NULL && PIPE[0] == '\0')
              fprintf(stderr,"Adding scaffolds from stdio ...\n");
            else
              fprintf(stderr,"Adding '%s.fasta' ...\n",core);
            fflush(stderr);
          }
    
        if (!append)
          flist[ofiles++] = core;

        // Check that the first line is a header line

        if (read[strlen(read)-1] != '\n')
          { fprintf(stderr,"File %s.fasta, Line 1: Fasta line is too long (> %d chars)\n",
                           core,MAX_NAME-2);
            goto error;
          }
        if (!eof && read[0] != '>')
          { fprintf(stderr,"File %s.fasta, Line 1: First header in fasta file is missing\n",core);
            goto error;
          }

        //  Read in all the sequences until end-of-file

        { int i, x, n;

          while (!eof)
            { int hlen;

              read[rlen] = '>';
              hlen = strlen(read+rlen);
              fwrite(read+rlen,1,hlen,hdrs);

              rlen  = 0;
              while (1)
                { eof = (fgets(read+rlen,MAX_NAME,input) == NULL);
                  nline += 1;
                  x = strlen(read+rlen)-1;
                  if (read[rlen+x] != '\n')
                    { fprintf(stderr,"File %s.fasta, Line %d:",core,nline);
                      fprintf(stderr," Fasta line is too long (> %d chars)\n",MAX_NAME-2);
                      goto error;
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

              n = 0;
              i = -1;
              while (i < rlen)
                { int pbeg, plen, clen;

                  while (i < rlen)
                    if (number[(int) read[++i]] < 4)
                      break;

                  if (i >= rlen) break;

                  pbeg = i;
                  prec.fpulse = pbeg;
                  prec.origin = n++;
                  prec.boff   = offset;
                  prec.coff   = hdrset;
                  prec.flags  = DB_BEST;
                  while (i < rlen)
                    { x = number[(int) read[i]];
                      if (x >= 4) break;
                      count[x] += 1;
                      read[i++] = (char) x;
                    }
                  prec.rlen = plen = i-pbeg;
                  ureads += 1;
                  totlen += plen;
                  if (plen > maxlen)
                    maxlen = plen;

                  Compress_Read(plen,read+pbeg);
                  clen = COMPRESSED_LEN(plen);
                  fwrite(read+pbeg,1,clen,bases);
                  offset += clen;

                  fwrite(&prec,sizeof(DAZZ_READ),1,indx);
                }
              hdrset += hlen;
            }
        }

        fprintf(ostub,DB_FDATA,ureads,core,core);

        if (PIPE == NULL)
          fclose(input);
        else
          break;
      }

    //  Update relevant fields in db record

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
        { fprintf(stderr,"%s: %s.dam is corrupted, read failed 2\n",Prog_Name,root);
          goto error;
        }
      dbpos = ftello(ostub);
      fprintf(ostub,DB_NBLOCK,0);
      if (fscanf(istub,DB_PARAMS,&size,&cutoff,&allflag) != 3)
        { fprintf(stderr,"%s: %s.dam is corrupted, read failed 3\n",Prog_Name,root);
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
            { fprintf(stderr,"%s: %s.dam is corrupted, read failed 4\n",Prog_Name,root);
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
          if (rlen >= cutoff)
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

  rewind(ostub);
  fprintf(ostub,DB_NFILE,ofiles);

  rewind(indx);
  fwrite(&db,sizeof(DAZZ_DB),1,indx);   //  Write the finalized db record into .idx

  if (istub != NULL)
    fclose(istub);
  fclose(ostub);
  fclose(indx);
  fclose(bases);
  fclose(hdrs);

  rename(Catenate(pwd,"/",root,".dbx"),dbname);   //  New image replaces old image

  exit (0);

  //  Error exit:  Either truncate or remove the .idx, .bps, and .hdr files as appropriate.
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
  if (hoff != 0)
    { fseeko(hdrs,0,SEEK_SET);
      if (ftruncate(fileno(hdrs),hoff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.hdr after error, truncate failed\n",
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
  if (hdrs != NULL)
    { fclose(hdrs);
      if (hoff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".hdr"));
    }
  if (ostub != NULL)
    { fclose(ostub);
      unlink(Catenate(pwd,"/",root,".dbx"));
    }
  if (istub != NULL)
    fclose(istub);

  exit (1);
}
