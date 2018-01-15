/********************************************************************************************
 *
 *  Recreate all the .fasta files that are in a specified DAM.
 *
 *  Author:  Gene Myers
 *  Date  :  May 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-vU] [-w<int(80)>] <path:dam>";

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  FILE       *dbfile, *hdrs;
  char       *dbfile_name, *hdrs_name;
  int         nfiles;
  int         VERBOSE, UPPER, WIDTH;

  //  Process arguments

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("DAM2fasta")

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

    UPPER   = 1 + flags['U'];
    VERBOSE = flags['v'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        exit (1);
      }
  }

  //  Open db

  { int   status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 0)
      { fprintf(stderr,"%s: Cannot be called on a .db: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (db->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
  }

  { char *pwd, *root;

    pwd    = PathTo(argv[1]);
    root   = Root(argv[1],".dam");
    dbfile_name = Strdup(Catenate(pwd,"/",root,".dam"),"Allocating db file name");
    hdrs_name   = Strdup(Catenate(pwd,PATHSEP,root,".hdr"),"Allocating header file name");
    dbfile = Fopen(dbfile_name,"r");
    hdrs   = Fopen(hdrs_name,"r");
    free(pwd);
    free(root);
    if (dbfile_name == NULL || hdrs_name == NULL || dbfile == NULL || hdrs == NULL)
      exit (1);
  }

  //  nfiles = # of files in data base

  FSCANF(dbfile,DB_NFILE,&nfiles)

  //  For each file do:

  { DAZZ_READ  *reads;
    char       *read;
    int         f, first;
    char        nstring[WIDTH+1];

    if (UPPER == 2)
      for (f = 0; f < WIDTH; f++)
        nstring[f] = 'N';
    else
      for (f = 0; f < WIDTH; f++)
        nstring[f] = 'n';
    nstring[WIDTH] = '\0';

    reads = db->reads;
    read  = New_Read_Buffer(db);
    first = 0;
    for (f = 0; f < nfiles; f++)
      { int   i, last, wpos;
        FILE *ofile;
        char  prolog[MAX_NAME], fname[MAX_NAME], header[MAX_NAME];

        //  Scan db image file line, create .fasta file for writing

        FSCANF(dbfile,DB_FDATA,&last,fname,prolog)

        if (strcmp(fname,"stdout") == 0)
          { ofile = stdout;

            if (VERBOSE)
              { fprintf(stderr,"Sending %d contigs to stdout ...\n",last-first);
                fflush(stdout);
              }
          }
        else
          { if ((ofile = Fopen(Catenate(".","/",fname,".fasta"),"w")) == NULL)
              exit (1);

            if (VERBOSE)
              { fprintf(stderr,"Creating %s.fasta ...\n",fname);
                fflush(stdout);
              }
          }

        //   For the relevant range of reads, write each to the file
        //     recreating the original headers with the index meta-data about each read

        wpos = 0;
        for (i = first; i < last; i++)
          { int        j, len, nlen, w;
            DAZZ_READ *r;

            r     = reads + i;
            len   = r->rlen;

            if (r->origin == 0)
              { if (i != first && wpos != 0)
                  { fprintf(ofile,"\n");
                    wpos = 0;
                  }
                FSEEKO(hdrs,r->coff,SEEK_SET)
                FGETS(header,MAX_NAME,hdrs)
                FPUTS(header,ofile)
              }

            if (r->fpulse != 0)
              { if (r->origin != 0)
                  nlen = r->fpulse - (reads[i-1].fpulse + reads[i-1].rlen);
                else
                  nlen = r->fpulse;

                for (j = 0; j+(w = WIDTH-wpos) <= nlen; j += w)
                  { FPRINTF(ofile,"%.*s\n",w,nstring)
                    wpos = 0;
                  }
                if (j < nlen)
                  { FPRINTF(ofile,"%.*s",nlen-j,nstring)
                    if (j == 0)
                      wpos += nlen;
                    else
                      wpos = nlen-j;
                  }
              }

            Load_Read(db,i,read,UPPER);

            for (j = 0; j+(w = WIDTH-wpos) <= len; j += w)
              { FPRINTF(ofile,"%.*s\n",w,read+j)
                wpos = 0;
              }
            if (j < len)
              { FPRINTF(ofile,"%s",read+j)
                if (j == 0)
                  wpos += len;
                else
                  wpos = len-j;
              }
          }
        if (wpos > 0)
          FPRINTF(ofile,"\n")
        if (ofile != stdout)
          FCLOSE(ofile)
        first = last;
      }
  }

  fclose(hdrs);
  fclose(dbfile);
  Close_DB(db);

  exit (0);
}
