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

static char *Usage = "[-af] [-x<int>] <path:db|dam>";

int main(int argc, char *argv[])
{ HITS_DB    db, dbs;
  int64      dbpos;
  FILE      *dbfile, *ixfile;
  int        nblocks;
  int        status;

  int        FORCE;
  int        ALL;
  int        CUTOFF;

  { int   i, j, k;
    int   flags[128];
    char *eptr;
    float size;

    ARG_INIT("DBtrim")

    CUTOFF = 0;
    size   = 200;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("af")
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

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
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

  { char    *pwd, *root;
    char     buffer[2*MAX_NAME+100];
    int      nfiles;
    int      all, cutoff;
    int64    size;
    int      i;

    pwd  = PathTo(argv[1]);
    if (status)
      { root   = Root(argv[1],".dam");
        dbfile = Fopen(Catenate(pwd,"/",root,".dam"),"r+");
      }
    else
      { root   = Root(argv[1],".db");
        dbfile = Fopen(Catenate(pwd,"/",root,".db"),"r+");
      }
    ixfile = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r+");
    if (dbfile == NULL || ixfile == NULL)
      exit (1);
    free(pwd);
    free(root);

    if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
      SYSTEM_ERROR
    for (i = 0; i < nfiles; i++)
      if (fgets(buffer,2*MAX_NAME+100,dbfile) == NULL)
        SYSTEM_ERROR

    if (fread(&dbs,sizeof(HITS_DB),1,ixfile) != 1)
      SYSTEM_ERROR

    if (dbs.cutoff >= 0)
      { if (!FORCE)
          { printf("You are about to reset the thresholds for the trimmed DB.\n");
            printf("This will invalidate any .las files produced by daligner\n");
            printf("Are you sure you want to proceed? [Y/N] ");
            fflush(stdout);
            if (fgets(buffer,100,stdin) == NULL)
              SYSTEM_ERROR
            if (index(buffer,'n') != NULL || index(buffer,'N') != NULL)
              { printf("Aborted\n");
                fflush(stdout);
                fclose(dbfile);
                exit (1);
              }
          }
      }
    else
      { fprintf(stderr,"%s: DB has not yet been split, use DBsplit\n",Prog_Name);
        exit (1);
      }

    fscanf(dbfile,DB_NBLOCK,&nblocks);

    dbpos = ftello(dbfile);
    fscanf(dbfile,DB_PARAMS,&size,&cutoff,&all);
    fseeko(dbfile,dbpos,SEEK_SET);
    fprintf(dbfile,DB_PARAMS,size,CUTOFF,ALL);
  }

  { HITS_READ *reads  = db.reads;
    int        uread, tread;
    int        rlen;
    int        b, u, t;

    u = 0;
    t = 0;
    fprintf(dbfile,DB_BDATA,0,0);
    for (b = 0; b < nblocks; b++)
      { dbpos = ftello(dbfile);
        fscanf(dbfile,DB_BDATA,&uread,&tread);

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

        fseeko(dbfile,dbpos,SEEK_SET);
        fprintf(dbfile,DB_BDATA,uread,t);
      }

    dbs.cutoff = CUTOFF;
    if (ALL)
      dbs.allarr |= DB_ALL;
    dbs.treads = t;
    rewind(ixfile);
    fwrite(&dbs,sizeof(HITS_DB),1,ixfile);
  }

  fclose(ixfile);
  fclose(dbfile);
  Close_DB(&db);

  exit (0);
}
