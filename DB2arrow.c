/********************************************************************************************
 *
 *  Recreate all the .arrow files that have been loaded into a specified database.
 *
 *  Author:  Gene Myers
 *  Date  :  October 2016
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "DB.h"

static char *Usage = "[-v] [-w<int(80)>] <path:db>";

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  FILE       *dbfile;
  char       *dbfile_name;
  int         VERBOSE, WIDTH;

  //  Process arguments

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("DB2arrow")

    WIDTH = 80;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vU")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        exit (1);
      }
  }

  //  Open db, and db stub file

  { int   status;
    char *pwd, *root;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 1)
      { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (db->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if ((db->allarr & DB_ARROW) == 0)
      { fprintf(stderr,"%s: There is no Arrow information in the DB: %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    pwd         = PathTo(argv[1]);
    root        = Root(argv[1],".db");
    dbfile_name = Strdup(Catenate(pwd,"/",root,".db"),"Allocting db file name");
    dbfile      = Fopen(dbfile_name,"r");
    free(pwd);
    free(root);
    if (dbfile_name == NULL || dbfile == NULL)
      exit (1);
  }

  //  For each cell do:

  { DAZZ_READ  *reads;
    char        lname[MAX_NAME];
    FILE       *ofile = NULL;
    int         f, first, last, ofirst, nfiles;
    char       *read;

    FSCANF(dbfile,DB_NFILE,&nfiles)

    reads = db->reads;
    read  = New_Read_Buffer(db);
    first = ofirst = 0;
    for (f = 0; f < nfiles; f++)
      { int   i;
        char  prolog[MAX_NAME], fname[MAX_NAME];

        //  Scan db image file line, create .arrow file for writing

        FSCANF(dbfile,DB_FDATA,&last,fname,prolog)

        if (f == 0 || strcmp(fname,lname) != 0)
          { if (f > 0)
              { if (ofile == stdout)
                  { fprintf(stderr," %d reads\n",first-ofirst);
                    fflush(stderr);
                  }
                else
                  FCLOSE(ofile)
              }

            if (strcmp(fname,"stdout") == 0)
              { ofile  = stdout;
                ofirst = first;

                if (VERBOSE)
                  { fprintf(stderr,"Sending to stdout ...");
                    fflush(stdout);
                  }
              }
            else
              { if ((ofile = Fopen(Catenate(".","/",fname,".arrow"),"w")) == NULL)
                  exit (1);

                if (VERBOSE)
                  { fprintf(stderr,"Creating %s.arrow ...\n",fname);
                    fflush(stdout);
                  }
              }

            strcpy(lname,fname);
          }

        //   For the relevant range of reads, write each to the file
        //     recreating the original headers with the index meta-data about each read

        for (i = first; i < last; i++)
          { int        j, len;
            uint64     big;
            float      snr[4];
            DAZZ_READ *r;

            r     = reads + i;
            len   = r->rlen;
            big   = *((uint64 *) &(r->coff));
            for (j = 0; j < 4; j++)
              { snr[3-j] = (big & 0xffff) / 100.;
                big    >>= 16;
              }

            FPRINTF(ofile,">%s",prolog)
            FPRINTF(ofile," SN=%.2f,%.2f,%.2f,%.2f",snr[0],snr[1],snr[2],snr[3])
            FPRINTF(ofile,"\n")

            Load_Arrow(db,i,read,1);

            for (j = 0; j+WIDTH < len; j += WIDTH)
              FPRINTF(ofile,"%.*s\n",WIDTH,read+j)
            if (j < len)
              FPRINTF(ofile,"%s\n",read+j)
          }

        first = last;
      }

    if (f > 0)
      { if (ofile == stdout)
          { fprintf(stderr," %d reads\n",first-ofirst);
            fflush(stderr);
          }
        else
          FCLOSE(ofile)
      }
  }

  fclose(dbfile);
  Close_DB(db);

  exit (0);
}
