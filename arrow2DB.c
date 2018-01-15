/*******************************************************************************************
 *
 *  Adds the given .arrow files to an existing DB "path".  The input files must be added in
 *  the same order as the .fasta files were and have the same root names, e.g. FOO.fasta
 *  and FOO.arrow.  The files can be added incrementally but must be added in the same order  
 *  as the .fasta files.  This is enforced by the program.  With the -l option set the
 *  compression scheme is a bit lossy to get more compression (see the description of dexqv
 *  in the DEXTRACTOR module).
 *
 *  Author:  Gene Myers
 *  Date  :  July 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "DB.h"
#include "QV.h"

//  Compiled in INTERACTIVE mode as all routines must return with an error
//    so that cleanup and restore is possible.

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-v] <path:db> ( -f<file> | -i | <input:arrow> ... )";

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
{ FILE      *istub;
  char      *root, *pwd;

  FILE      *arrow, *indx;
  int64      boff;

  DAZZ_DB    db;
  DAZZ_READ *reads;
  int        nfiles;

  int        VERBOSE;
  int        PIPE;
  FILE      *INFILE;

  //  Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("arrow2DB")

    INFILE = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vli")
            break;
          case 'f':
            INFILE = fopen(argv[i]+2,"r");
            if (INFILE == NULL)
              { fprintf(stderr,"%s: Cannot open file of inputs '%s'\n",Prog_Name,argv[i]+2);
                exit (1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    PIPE    = flags['i'];

    if (INFILE != NULL && PIPE)
      { fprintf(stderr,"%s: Cannot use both -f and -i together\n",Prog_Name);
        exit (1);
      }

    if ( (INFILE == NULL && ! PIPE && argc <= 2) || 
        ((INFILE != NULL || PIPE) && argc != 2))
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"      -f: import files listed 1/line in given file.\n");
        fprintf(stderr,"      -i: import data from stdin.\n");
        fprintf(stderr,"        : otherwise, import sequence of specified files.\n");
        exit (1);
      }
  }

  //  Open DB stub file, index, and .arw file for appending.  Load db and read records,
  //    get number of cells from stub file, and note current offset to end of .arw

  root   = Root(argv[1],".db");
  pwd    = PathTo(argv[1]);
  istub  = Fopen(Catenate(pwd,"/",root,".db"),"r");
  if (istub == NULL)
    exit (1);
  if (fscanf(istub,DB_NFILE,&nfiles) != 1)
    { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
      exit (1);
    }

  indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r+");
  if (indx == NULL)
    exit (1);
  if (fread(&db,sizeof(DAZZ_DB),1,indx) != 1)
    { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
      exit (1);
    }

  reads = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*db.ureads,"Allocating DB index");
  if (reads == NULL)
    exit (1);
  if (fread(reads,sizeof(DAZZ_READ),db.ureads,indx) != (size_t) (db.ureads))
    { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
      exit (1);
    }

  if (reads[0].coff >= 0 && (db.allarr & DB_ARROW) == 0)
    { fprintf(stderr,"%s: Database %s has Quiver data!\n",Prog_Name,root);
      exit (1);
    }

  arrow = NULL;
  boff  = 0;

  if (reads[0].coff < 0)
    arrow = Fopen(Catenate(pwd,PATHSEP,root,".arw"),"w");
  else
    arrow = Fopen(Catenate(pwd,PATHSEP,root,".arw"),"r+");

  if (arrow == NULL)
    goto error;
  fseeko(arrow,0,SEEK_END);
  boff = ftello(arrow);

  //  Do a merged traversal of cell lines in .db stub file and .arrow files to be
  //    imported, driving the loop with the cell line #

  { FILE          *input = NULL;
    char          *path = NULL;
    char          *core = NULL;
    char          *read;
    int            rmax, rlen, eof;
    File_Iterator *ng = NULL;
    char           lname[MAX_NAME];
    int            first, last, cline;
    int            cell;

    //  Buffer for accumulating .arrow sequence over multiple lines

    rmax  = MAX_NAME + 60000;
    read  = (char *) Malloc(rmax+1,"Allocating line buffer");
    if (read == NULL)
      goto error;

    if (!PIPE)
      { ng = init_file_iterator(argc,argv,INFILE,2);
        if (ng == NULL)
          goto error;
      }

    eof = 0;
    for (cell = 0; cell < nfiles; cell++)
      { char  prolog[MAX_NAME], fname[MAX_NAME];

        if (cell == 0)

          //  First addition, a pipe: find the first cell that does not have .arrow's yet
          //     (error if none) and set input source to stdin.

          if (PIPE)
            { first = 0;
              while (cell < nfiles)
                { if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
                    { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",core,Prog_Name);
                      goto error;
                    }
                  if (reads[first].coff < 0)
                    break;
                  first = last;
                  cell += 1;
                }
              if (cell >= nfiles)
                { fprintf(stderr,"%s: All .arrows's have already been added !?\n",Prog_Name);
                  goto error;
                }

              input = stdin;

              if (VERBOSE)
                { fprintf(stderr,"Adding arrows's from stdin ...\n");
                  fflush(stderr);
                }
              cline = 0;
            }

          //  First addition, not a pipe: then get first .arrow file name (error if not one) to
          //    add, find the first cell name whose file name matches (error if none), check that
          //    the previous .arrow's have been added and this is the next slot.  Then open
          //    the .arrow file for compression

          else
            { if (! next_file(ng))
                { fprintf(stderr,"%s: file list is empty!\n",Prog_Name);
                  goto error;
                }
              if (ng->name == NULL)
                goto error;
  
              core = Root(ng->name,".arrow");
              path = PathTo(ng->name);
              if ((input = Fopen(Catenate(path,"/",core,".arrow"),"r")) == NULL)
                goto error;
  
              first = 0;
              while (cell < nfiles)
                { if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
                    { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",core,Prog_Name);
                      goto error;
                    }
                  if (strcmp(core,fname) == 0)
                    break;
                  first = last;
                  cell += 1;
                }
              if (cell >= nfiles)
                { fprintf(stderr,"%s: %s.fasta has never been added to DB\n",Prog_Name,core);
                  goto error;
                }
        
              if (first > 0 && reads[first-1].coff < 0)
                { fprintf(stderr,"%s: Predecessor of %s.arrow has not been added yet\n",
                                 Prog_Name,core);
                  goto error;
                }
              if (reads[first].coff >= 0)
                { fprintf(stderr,"%s: %s.arrow has already been added\n",Prog_Name,core);
                  goto error;
                }
  
              if (VERBOSE)
                { fprintf(stderr,"Adding '%s.arrow' ...\n",core);
                  fflush(stderr);
                }
              cline = 0;
            }

        //  Not the first addition: get next cell line.  If not a pipe and the file name is new,
        //    then close the current .arrow, open the next one and after ensuring the names
        //    match, open it for incorporation

        else
          { first = last;  
            strcpy(lname,fname);
            if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
              { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",core,Prog_Name);
                goto error;
              }
            if (PIPE)
              { int c;
                if ((c = fgetc(input)) == EOF)
                  break;
                ungetc(c,input);
              }
            else if (strcmp(lname,fname) != 0)
              { if ( ! eof)
                  { fprintf(stderr,"%s: Too many reads in %s.arrow while handling %s.fasta\n",
                                   Prog_Name,core,fname);
                    goto error;
                  }

                fclose(input);
                free(path);
                free(core);

                if ( ! next_file(ng))
                  break;
                if (ng->name == NULL)
                  goto error;

                path = PathTo(ng->name);
                core = Root(ng->name,".arrow");
                if ((input = Fopen(Catenate(path,"/",core,".arrow"),"r")) == NULL)
                  goto error;

                if (strcmp(core,fname) != 0)
                  { fprintf(stderr,"%s: Files not being added in order (expect %s, given %s)\n",
                                   Prog_Name,fname,core);
                    goto error;
                  }

                if (VERBOSE)
                  { fprintf(stderr,"Adding '%s.arrow' ...\n",core);
                    fflush(stderr);
                  }
                cline = 0;
              }
          }

        //  If first cell or source is a new file, then start IO

        if (cline == 0)
          {
            // Read in first line and make sure it is a header in PACBIO format.

            rlen = 0;
            eof  = (fgets(read,MAX_NAME,input) == NULL);
            if (read[strlen(read)-1] != '\n')
              { fprintf(stderr,"File %s.arrow, Line 1: Fasta line is too long (> %d chars)\n",
                               core,MAX_NAME-2);
                goto error;
              }
            if (!eof && read[0] != '>')
              { fprintf(stderr,"File %s.arrow, Line 1: First header in arrow file is missing\n",
                                core);
                goto error;
              }
          }

        //  Compress reads [first..last) from open .arrow appending to .arw and record
        //    snr in .coff field of reads (offset is the same as for the DNA sequence, .boff)

        { int i, x;

          for (i = first; i < last; i++)
            { char  *find;
              int    clen;
              float  snr[4];
              uint16 cnr[4];

              if (eof)
                { if (PIPE)
                    fprintf(stderr,"%s: Insufficient # of reads on input while handling %s.arrow\n",
                                   Prog_Name,fname);
                  else
                    { fprintf(stderr,"%s: Insufficient # of reads in %s.arrow while handling",
                                     Prog_Name,core);
                      fprintf(stderr," %s.arrow\n",fname);
                    }
                  goto error;
                }

              find = index(read+(rlen+1),' ');
              if (find == NULL)
                { fprintf(stderr,"File %s.arrow, Line %d: Pacbio header line format error\n",
                                 core,cline);
                  goto error;
                }
              *find = '\0';
              if (strcmp(read+(rlen+1),prolog) != 0)
                { fprintf(stderr,"File %s.arrow, Line %d: Pacbio prolog doesn't match DB entry\n",
                                 core,cline);
                  goto error;
                }
              *find = ' ';
              x = sscanf(find+1," SN=%f,%f,%f,%f\n",snr,snr+1,snr+2,snr+3);
              if (x != 4)
                { fprintf(stderr,"File %s.arrow, Line %d: Pacbio header line format error\n",
                                 core,cline);
                  goto error;
                }

              rlen  = 0;
              while (1)
                { eof = (fgets(read+rlen,MAX_NAME,input) == NULL);
                  cline += 1;
                  x = strlen(read+rlen)-1;
                  if (read[rlen+x] != '\n')
                    { if (read[rlen] == '>')
                        { fprintf(stderr,"File %s.arrow, Line %d:",core,cline);
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
                    { rmax = ((int) (1.2 * rmax)) + 1000 + MAX_NAME;
                      read = (char *) realloc(read,rmax+1);
                      if (read == NULL)
                        { fprintf(stderr,"File %s.arrow, Line %d:",core,cline);
                          fprintf(stderr," Out of memory (Allocating line buffer)\n");
                          goto error;
                        }
                    }
                }
              read[rlen] = '\0';

              for (x = 0; x < 4; x++)
                cnr[x] = (uint32) (snr[x] * 100.);

              *((uint64  *) &(reads[i].coff)) = ((uint64) cnr[0]) << 48 |
                                                ((uint64) cnr[1]) << 32 |
                                                ((uint64) cnr[2]) << 16 |
                                                ((uint64) cnr[3]);

              Number_Arrow(read);
              Compress_Read(rlen,read);
              clen = COMPRESSED_LEN(rlen);
              fwrite(read,1,clen,arrow);
            }
        }
      }

    if (!eof)
      { if (PIPE)
          fprintf(stderr,"%s: Too many reads on input while handling %s.fasta\n",
                         Prog_Name,lname);
        else
          fprintf(stderr,"%s: Too many reads in %s.arrow while handling %s.fasta\n",
                         Prog_Name,core,lname);
        goto error;
      }
    if ( ! PIPE && cell >= nfiles)
      { fclose(input);
        free(core);
        free(path);
        if (next_file(ng))
          { if (ng->name == NULL)
              goto error;
            core = Root(ng->name,".arrow");
            fprintf(stderr,"%s: %s.fasta has never been added to DB\n",Prog_Name,core);
            goto error;
          }
      }
  }

  //  Write the db record and read index into .idx and clean up

  db.allarr |= DB_ARROW;
  rewind(indx);
  fwrite(&db,sizeof(DAZZ_DB),1,indx);
  fwrite(reads,sizeof(DAZZ_READ),db.ureads,indx);

  fclose(istub);
  fclose(indx);
  fclose(arrow);

  exit (0);

  //  Error exit:  Either truncate or remove the .arw file as appropriate.

error:
  if (boff != 0)
    { fseeko(arrow,0,SEEK_SET);
      if (ftruncate(fileno(arrow),boff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.arw after error, truncate failed\n",
                       Prog_Name,root);
    }
  if (arrow != NULL)
    { fclose(arrow);
      if (boff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".arw"));
    }
  fclose(istub);
  fclose(indx);

  exit (1);
}
