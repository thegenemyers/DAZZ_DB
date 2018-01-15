/*******************************************************************************************
 *
 *  Adds the given .quiva files to an existing DB "path".  The input files must be added in
 *  the same order as the .fasta files were and have the same root names, e.g. FOO.fasta
 *  and FOO.quiva.  The files can be added incrementally but must be added in the same order  
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

static char *Usage = "[-v] <path:db> ( -f<file> | -i | <input:quiva> ... )";

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

  FILE      *quiva, *indx;
  int64      coff;

  DAZZ_DB    db;
  DAZZ_READ *reads;
  int        nfiles;

  FILE      *temp;
  char      *tname;

  int        VERBOSE;
  int        PIPE;
  FILE      *INFILE;

  //  Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("quiva2DB")

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
        fprintf(stderr,"\n");
        fprintf(stderr,"      -f: import files listed 1/line in given file.\n");
        fprintf(stderr,"      -i: import data from stdin.\n");
        fprintf(stderr,"        : otherwise, import sequence of specified files.\n");
        exit (1);
      }
  }

  //  Open DB stub file, index, and .qvs file for appending.  Load db and read records,
  //    get number of cells from stub file, and note current offset to end of .qvs

  root   = Root(argv[1],".db");
  pwd    = PathTo(argv[1]);
  istub  = Fopen(Catenate(pwd,"/",root,".db"),"r");
  if (istub == NULL)
    { fprintf(stderr,"%s",Ebuffer);
      exit (1);
    }
  if (fscanf(istub,DB_NFILE,&nfiles) != 1)
    { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
      exit (1);
    }

  indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r+");
  if (indx == NULL)
    { fprintf(stderr,"%s",Ebuffer);
      exit (1);
    }
  if (fread(&db,sizeof(DAZZ_DB),1,indx) != 1)
    { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
      exit (1);
    }
  if ((db.allarr & DB_ARROW) != 0)
    { fprintf(stderr,"%s: Database %s has Arrow data!\n",Prog_Name,root);
      exit (1);
    }

  reads = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*db.ureads,"Allocating DB index");
  if (reads == NULL)
    { fprintf(stderr,"%s",Ebuffer);
      exit (1);
    }
  if (fread(reads,sizeof(DAZZ_READ),db.ureads,indx) != (size_t) (db.ureads))
    { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
      exit (1);
    }

  quiva = NULL;
  temp  = NULL;
  coff  = 0;

  if (reads[0].coff < 0)
    quiva = Fopen(Catenate(pwd,PATHSEP,root,".qvs"),"w");
  else
    quiva = Fopen(Catenate(pwd,PATHSEP,root,".qvs"),"r+");

  tname = Strdup(Catenate(".",PATHSEP,root,Numbered_Suffix("",getpid(),".tmp")),
                 "Allocating temporary name");
  temp = Fopen(tname,"w+");

  if (quiva == NULL || temp == NULL)
    { fprintf(stderr,"%s",Ebuffer);
      goto error;
    }
  fseeko(quiva,0,SEEK_END);
  coff = ftello(quiva);

  //  Do a merged traversal of cell lines in .db stub file and .quiva files to be
  //    imported, driving the loop with the cell line #

  { FILE          *input = NULL;
    char          *path = NULL;
    char          *core = NULL;
    File_Iterator *ng = NULL;
    char           lname[MAX_NAME];
    int            first, last, cline;
    int            cell;

    if (!PIPE)
      { ng = init_file_iterator(argc,argv,INFILE,2);
        if (ng == NULL)
          { fprintf(stderr,"%s",Ebuffer);
            goto error;
          }
      }

    for (cell = 0; cell < nfiles; cell++)
      { char  prolog[MAX_NAME], fname[MAX_NAME];

        if (cell == 0)

          //  First addition, a pipe: find the first cell that does not have .quiva's yet
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
                { fprintf(stderr,"%s: All .quiva's have already been added !?\n",Prog_Name);
                  goto error;
                }

              input = stdin;

              if (VERBOSE)
                { fprintf(stderr,"Adding quiva's from stdin ...\n");
                  fflush(stderr);
                }
              cline = 0;
            }

          //  First addition, not a pipe: then get first .quiva file name (error if not one) to
          //    add, find the first cell name whose file name matches (error if none), check that
          //    the previous .quiva's have been added and this is the next slot.  Then open
          //    the .quiva file for compression

          else
            { if (! next_file(ng))
                { fprintf(stderr,"%s: file list is empty!\n",Prog_Name);
                  goto error;
                }
              if (ng->name == NULL)
                { fprintf(stderr,"%s",Ebuffer);
                  goto error;
                }
  
              core = Root(ng->name,".quiva");
              path = PathTo(ng->name);
              if ((input = Fopen(Catenate(path,"/",core,".quiva"),"r")) == NULL)
                { fprintf(stderr,"%s",Ebuffer);
                  goto error;
                }
  
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
                { fprintf(stderr,"%s: Predecessor of %s.quiva has not been added yet\n",
                                 Prog_Name,core);
                  goto error;
                }
              if (reads[first].coff >= 0)
                { fprintf(stderr,"%s: %s.quiva has already been added\n",Prog_Name,core);
                  goto error;
                }
  
              if (VERBOSE)
                { fprintf(stderr,"Adding '%s.quiva' ...\n",core);
                  fflush(stderr);
                }
              cline = 0;
            }

        //  Not the first addition: get next cell line.  If not a pipe and the file name is new,
        //    then close the current .quiva, open the next one and after ensuring the names
        //    match, open it for compression

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
              { if (fgetc(input) != EOF)
                  { fprintf(stderr,"%s: Too many reads in %s.quiva while handling %s.fasta\n",
                                   Prog_Name,core,fname);
                    goto error;
                  }

                fclose(input);
                free(path);
                free(core);

                if ( ! next_file(ng))
                  break;
                if (ng->name == NULL)
                  { fprintf(stderr,"%s",Ebuffer);
                    goto error;
                  }

                path = PathTo(ng->name);
                core = Root(ng->name,".quiva");
                if ((input = Fopen(Catenate(path,"/",core,".quiva"),"r")) == NULL)
                  { fprintf(stderr,"%s",Ebuffer);
                    goto error;
                  }

                if (strcmp(core,fname) != 0)
                  { fprintf(stderr,"%s: Files not being added in order (expect %s, given %s)\n",
                                   Prog_Name,fname,core);
                    goto error;
                  }

                if (VERBOSE)
                  { fprintf(stderr,"Adding '%s.quiva' ...\n",core);
                    fflush(stderr);
                  }
                cline = 0;
              }
          }

        //  Compress reads [first..last) from open .quiva appending to .qvs and record
        //    offset in .coff field of reads (offset of first in a cell is to the compression
        //    table).

        { int64     qpos;
          QVcoding *coding;
          int       i, s;

          rewind(temp);
          if (ftruncate(fileno(temp),0) < 0)
            { fprintf(stderr,"%s: System error: could not truncate temporary file\n",Prog_Name);
              goto error;
            }
          Set_QV_Line(cline);
          s = QVcoding_Scan(input,last-first,temp);
          if (s < 0)
            { fprintf(stderr,"%s",Ebuffer);
              goto error;
            }
          if (s != last-first)
            { if (PIPE)
                fprintf(stderr,"%s: Insufficient # of reads on input while handling %s.fasta\n",
                               Prog_Name,fname);
              else
                fprintf(stderr,"%s: Insufficient # of reads in %s.quiva while handling %s.fasta\n",
                               Prog_Name,core,fname);
              goto error;
            }

          coding = Create_QVcoding(0);
          if (coding == NULL)
            { fprintf(stderr,"%s",Ebuffer);
              goto error;
            }

          coding->prefix = Strdup(".qvs","Allocating header prefix");
          if (coding->prefix == NULL)
            { fprintf(stderr,"%s",Ebuffer);
              goto error;
            }

          qpos = ftello(quiva);
          Write_QVcoding(quiva,coding);

          //  Then compress and append to the .qvs each compressed QV entry
 
          rewind(temp);
          Set_QV_Line(cline);
          for (i = first; i < last; i++)
            { s = Read_Lines(temp,1);
              if (s < -1)
                { fprintf(stderr,"%s",Ebuffer);
                  goto error;
                }
              reads[i].coff = qpos;
              s = Compress_Next_QVentry(temp,quiva,coding,0);
              if (s < 0)
                { fprintf(stderr,"%s",Ebuffer);
                  goto error;
                }
              if (s != reads[i].rlen)
                { fprintf(stderr,"%s: Length of quiva %d is different than fasta in DB\n",
                                 Prog_Name,i+1);
                  goto error;
                }
              qpos = ftello(quiva);
            }
          cline = Get_QV_Line();

          Free_QVcoding(coding);
        }
      }

    if (fgetc(input) != EOF)
      { if (PIPE)
          fprintf(stderr,"%s: Too many reads on input while handling %s.fasta\n",
                         Prog_Name,lname);
        else
          fprintf(stderr,"%s: Too many reads in %s.quiva while handling %s.fasta\n",
                         Prog_Name,core,lname);
        goto error;
      }
    if ( ! PIPE && cell >= nfiles)
      { fclose(input);
        free(core);
        free(path);
        if (next_file(ng))
          { if (ng->name == NULL)
              { fprintf(stderr,"%s",Ebuffer);
                goto error;
              }
            core = Root(ng->name,".quiva");
            fprintf(stderr,"%s: %s.fasta has never been added to DB\n",Prog_Name,core);
            goto error;
          }
      }
  }

  //  Write the db record and read index into .idx and clean up

  rewind(indx);
  fwrite(&db,sizeof(DAZZ_DB),1,indx);
  fwrite(reads,sizeof(DAZZ_READ),db.ureads,indx);

  fclose(istub);
  fclose(indx);
  fclose(quiva);
  fclose(temp);
  unlink(tname);

  exit (0);

  //  Error exit:  Either truncate or remove the .qvs file as appropriate.

error:
  if (coff != 0)
    { fseeko(quiva,0,SEEK_SET);
      if (ftruncate(fileno(quiva),coff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.qvs after error, truncate failed\n",
                       Prog_Name,root);
    }
  if (quiva != NULL)
    { fclose(quiva);
      if (coff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".qvs"));
    }
  if (temp != NULL)
    { fclose(temp);
      unlink(tname);
    }
  fclose(istub);
  fclose(indx);

  exit (1);
}
