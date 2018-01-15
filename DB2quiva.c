/********************************************************************************************
 *
 *  Recreate all the .quiva files that have been loaded into a specified database.
 *
 *  Author:  Gene Myers
 *  Date  :  May 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "DB.h"
#include "QV.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-vU] <path:db>";

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  FILE       *dbfile, *quiva;
  char       *dbfile_name;
  int         VERBOSE, UPPER;

  //  Process arguments

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("DB2quiva")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vU") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    UPPER   = flags['U'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        exit (1);
      }
  }

  //  Open db, db stub file, and .qvs file

  { char *pwd, *root;
    int   status;

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
    if (db->reads[0].coff < 0 || (db->allarr & DB_ARROW) != 0)
      { fprintf(stderr,"%s: There is no Quiver information in the DB: %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    pwd         = PathTo(argv[1]);
    root        = Root(argv[1],".db");
    dbfile_name = Strdup(Catenate(pwd,"/",root,".db"),"Allocating db file name");
    dbfile      = Fopen(dbfile_name,"r");
    quiva       = Fopen(Catenate(pwd,PATHSEP,root,".qvs"),"r");
    free(pwd);
    free(root);
    if (dbfile_name == NULL || dbfile == NULL || quiva == NULL)
      exit (1);
  }

  //  For each cell do:

  { DAZZ_READ  *reads;
    char        lname[MAX_NAME];
    FILE       *ofile = NULL;
    int         f, first, last, ofirst, nfiles;
    QVcoding   *coding;
    char      **entry;

    FSCANF(dbfile,DB_NFILE,&nfiles)

    reads = db->reads;
    entry = New_QV_Buffer(db);
    first = ofirst = 0;
    for (f = 0; f < nfiles; f++)
      { int   i;
        char  prolog[MAX_NAME], fname[MAX_NAME];

        //  Scan db image file line, create .quiva file for writing

        if (reads[first].coff < 0) break;

        FSCANF(dbfile,DB_FDATA,&last,fname,prolog)

        if (f == 0 || strcmp(fname,lname) != 0)
          { if (f > 0)
              { if (ofile == stdout)
                  { fprintf(stderr," %d quivas\n",first-ofirst);
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
              { if ((ofile = Fopen(Catenate(".","/",fname,".quiva"),"w")) == NULL)
                  exit (1);

                if (VERBOSE)
                  { fprintf(stderr,"Creating %s.quiva ...\n",fname);
                    fflush(stderr);
                  }
              }

            strcpy(lname,fname);
          }

        //   For the relevant range of reads, write the header for each to the file
        //     and then uncompress and write the quiva entry for each

        coding = Read_QVcoding(quiva);

        for (i = first; i < last; i++)
          { int        e, flags, qv, rlen;
            DAZZ_READ *r;

            r     = reads + i;
            flags = r->flags;
            rlen  = r->rlen;
            qv    = (flags & DB_QV);
            FPRINTF(ofile,"@%s/%d/%d_%d",prolog,r->origin,r->fpulse,r->fpulse+rlen)
            if (qv > 0)
              FPRINTF(ofile," RQ=0.%3d",qv)
            FPRINTF(ofile,"\n")

            Uncompress_Next_QVentry(quiva,entry,coding,rlen);

            if (UPPER)
              { char *deltag = entry[1];
                int   j;

                for (j = 0; j < rlen; j++)
                  deltag[j] -= 32;
              }

            for (e = 0; e < 5; e++)
              FPRINTF(ofile,"%.*s\n",rlen,entry[e])
          }

        first = last;
      }

    if (f > 0)
      { if (ofile == stdout)
          { fprintf(stderr," %d quivas\n",first-ofirst);
            fflush(stderr);
          }
        else
          FCLOSE(ofile)
      }
  }

  fclose(quiva);
  fclose(dbfile);
  Close_DB(db);

  exit (0);
}
