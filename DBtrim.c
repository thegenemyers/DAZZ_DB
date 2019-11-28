/*******************************************************************************************
 *
 *  Reset the trimming parameters for a .db:
 *     Rewrite the .db or .dam file with the new thresholds and the new read counts for
 *     each trimmed block.
 *
 *  Author:  Gene Myers
 *  Date  :  September 2017
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-aflm] [-x<int>] <path:db|dam>";

static DAZZ_READ *reads;

static int MSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (reads[x].rlen - reads[y].rlen);
}

int main(int argc, char *argv[])
{ DAZZ_DB    db, dbs;
  int64      dbpos;
  FILE      *dbfile, *ixfile;
  char      *dbfile_name, *ixfile_name;
  int        nblocks;
  int        status;
  int        nreads;

  int        medmax = 0;
  int       *msort = NULL;

  int        FORCE;
  int        ALL;
  int        CUTOFF;
  int        SELECT;

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("DBtrim")

    SELECT = -1;
    CUTOFF = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("aflm")
            break;
          case 'x':
            ARG_NON_NEGATIVE(CUTOFF,"Min read length cutoff")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    ALL   = flags['a'];
    FORCE = flags['f'];

    if (flags['l'])
      { SELECT = 0;
        if (flags['m'])
          { fprintf(stderr,"%s: Cannot set both -l and -m, mutually exclusive\n",Prog_Name);
            exit (1);
          }
      }
    else if (flags['m'])
      SELECT = 1;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -x: Trimmed DB has reads >= this threshold.\n");
        fprintf(stderr,"      -a: Trimmed DB contains all reads from a well (not just longest).\n");
        fprintf(stderr,"      -f: Force the new trim setting even if already set.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: Set primary read for a well to be the longest.\n");
        fprintf(stderr,"      -m: Set primary read for a well to be the median.\n");
        exit (1);
      }
  }

  //  Open db

  status = Open_DB(argv[1],&db);
  if (status < 0)
    exit (1);
  if (db.part > 0)
    { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  reads  = db.reads;
  nreads = db.ureads;

  { char    *pwd, *root;
    char     buffer[2*MAX_NAME+100];
    int      nfiles;
    int      all, cutoff;
    int64    size;
    int      i;

    pwd  = PathTo(argv[1]);
    if (status)
      { root   = Root(argv[1],".dam");
        dbfile_name = Strdup(Catenate(pwd,"/",root,".dam"),"Allocating db file name");
      }
    else
      { root   = Root(argv[1],".db");
        dbfile_name = Strdup(Catenate(pwd,"/",root,".db"),"Allocating db file name");
      }
    ixfile_name = Strdup(Catenate(pwd,PATHSEP,root,".idx"),"Allocating index file name");
    dbfile = Fopen(dbfile_name,"r+");
    ixfile = Fopen(ixfile_name,"r+");
    if (dbfile_name == NULL || ixfile_name == NULL || dbfile == NULL || ixfile == NULL)
      exit (1);
    free(pwd);
    free(root);

    FSCANF(dbfile,DB_NFILE,&nfiles)
    for (i = 0; i < nfiles; i++)
      FGETS(buffer,2*MAX_NAME+100,dbfile)

    FFREAD(&dbs,sizeof(DAZZ_DB),1,ixfile)

    if (dbs.cutoff >= 0)
      { if (!FORCE)
          { printf("You are about to reset the thresholds for the trimmed DB.\n");
            printf("This will invalidate any .las files produced by daligner\n");
            printf("Are you sure you want to proceed? [Y/N] ");
            fflush(stdout);
            fgets(buffer,100,stdin);
            if (index(buffer,'n') != NULL || index(buffer,'N') != NULL)
              { printf("Aborted\n");
                fflush(stdout);
                fclose(ixfile);
                fclose(dbfile);
                exit (1);
              }
          }
      }
    else
      { fprintf(stderr,"%s: DB has not yet been split, use DBsplit\n",Prog_Name);
        exit (1);
      }

    FSCANF(dbfile,DB_NBLOCK,&nblocks)

    FTELLO(dbpos,dbfile)
    FSCANF(dbfile,DB_PARAMS,&size,&cutoff,&all)
    FSEEKO(dbfile,dbpos,SEEK_SET)
    FPRINTF(dbfile,DB_PARAMS,size,CUTOFF,ALL)
  }

  if (SELECT >= 0)
    { int        i, j, k, h;
      int        off;

      off = (~ DB_BEST);

      for (i = 0; i < nreads; i = j)
        { j = i+1;
          reads[i].flags &= off;
          while ((reads[j].flags & DB_CCS) != 0)
            reads[j++].flags &= off;

          if (j-i <= 1)
            { reads[i].flags |= DB_BEST;
              continue;
            }

          if (SELECT == 0)
            { int mlen;

              mlen = reads[h = i].rlen;
              for (k = i+1; k < j; k++)
                if (mlen < reads[k].rlen)
                  mlen = reads[h = k].rlen;
            }
          else
            { if (j-i > medmax)
                { medmax = 1.2*(j-i) + 1000;
                  msort  = (int *) Realloc(msort,sizeof(int)*medmax,"Allocating Median Array");
                  if (msort == NULL)
                    exit (1);
                }
              for (k = i; k < j; k++)
                msort[k-i] = k;
              qsort(msort,j-i,sizeof(int),MSORT);
              h = msort[(j-i)/2];
            }
          reads[h].flags |= DB_BEST;
        }
    }

  { DAZZ_READ *reads  = db.reads;
    int        uread, tread;
    int        rlen;
    int        b, u, t;

    u = 0;
    t = 0;
    fprintf(dbfile,DB_BDATA,0,0);
    for (b = 0; b < nblocks; b++)
      { FTELLO(dbpos,dbfile);
        FSCANF(dbfile,DB_BDATA,&uread,&tread)

        if (ALL)
          while (u < uread)
            { rlen = reads[u++].rlen;
              if (rlen >= CUTOFF)
                t += 1;
            }
        else
          while (u < uread)
            { rlen = reads[u].rlen;
              if (rlen >= CUTOFF && (reads[u].flags & DB_BEST) != 0)
                t += 1;
              u += 1;
            }

        FSEEKO(dbfile,dbpos,SEEK_SET)
        FPRINTF(dbfile,DB_BDATA,uread,t)
      }

    dbs.cutoff = CUTOFF;
    if (ALL)
      dbs.allarr |= DB_ALL;
    else
      dbs.allarr &= ~DB_ALL;
    dbs.treads = t;
    FSEEKO(ixfile,0,SEEK_SET)
    FFWRITE(&dbs,sizeof(DAZZ_DB),1,ixfile)
    if (SELECT >= 0)
      FFWRITE(reads,sizeof(DAZZ_READ),nreads,ixfile)
  }

  FCLOSE(ixfile)
  FCLOSE(dbfile)
  Close_DB(&db);

  exit (0);
}
