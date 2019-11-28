/*******************************************************************************************
 *
 *  Split a .db into a set of sub-database blocks for use by the Dazzler:
 *     Divide the database <path>.db conceptually into a series of blocks referable to on the
 *     command line as <path>.1.db, <path>.2.db, ...  If the -x option is set then all reads
 *     less than the given length are ignored, and if the -a option is not set then secondary
 *     reads from a given well are also ignored.  The remaining reads are split amongst the
 *     blocks so that each block is of size -s * 1Mbp except for the last which necessarily
 *     contains a smaller residual.  The default value for -s is 400Mbp because blocks of this
 *     size can be compared by our "overlapper" dalign in roughly 16Gb of memory.  The blocks
 *     are very space efficient in that their sub-index of the master .idx is computed on the
 *     fly when loaded, and the .bps file of base pairs is shared with the master DB.  Any
 *     tracks associated with the DB are also computed on the fly when loading a database block.
 *
 *  Author:  Gene Myers
 *  Date  :  September 2013
 *  Mod   :  New splitting definition to support incrementality, and new stub file format
 *  Date  :  April 2014
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

static char *Usage = "[-aflm] [-x<int>] [-s<double(200.)>] <path:db|dam>";

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
  int        status;
  int        nreads;

  int        medmax = 0;
  int       *msort = NULL;

  int        FORCE;
  int        ALL;
  int        CUTOFF;
  int64      SIZE;
  int        SELECT;

  { int   i, j, k;
    int   flags[128];
    char *eptr;
    float size;

    ARG_INIT("DBsplit")

    SELECT = -1;
    CUTOFF = 0;
    size   = 200;

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
          case 's':
            ARG_REAL(size)
            if (size <= 0.)
              { fprintf(stderr,"%s: Block size must be a positive number\n",Prog_Name);
                exit (1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    SIZE  = size*1000000ll;
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
        fprintf(stderr,"      -s: Target size of blocks (in Mbp).\n");
        fprintf(stderr,"      -x: Trimmed DB has reads >= this threshold.\n");
        fprintf(stderr,"      -a: Trimmed DB contains all reads from a well (not just longest).\n");
        fprintf(stderr,"      -f: Force the split to occur even if already split.\n");
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
    int      i;

    pwd    = PathTo(argv[1]);
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

    if (dbs.cutoff >= 0 && !FORCE)
      { printf("You are about to overwrite the current partition settings.  This\n");
        printf("will invalidate any tracks, overlaps, and other derivative files.\n");
        printf("Are you sure you want to proceed? [Y/N] ");
        fflush(stdout);
        fgets(buffer,100,stdin);
        if (index(buffer,'n') != NULL || index(buffer,'N') != NULL)
          { printf("Aborted\n");
            fflush(stdout);
            fclose(dbfile);
            fclose(ixfile);
            exit (1);
          }
      }

    FTELLO(dbpos,dbfile)
    FSEEKO(dbfile,dbpos,SEEK_SET)
    FPRINTF(dbfile,DB_NBLOCK,0)
    FPRINTF(dbfile,DB_PARAMS,SIZE,CUTOFF,ALL)
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


  { int64      totlen;
    int        nblock, ireads, treads, rlen, fno;
    int        i, upos, css;

    css    = 0;
    nblock = 0;
    totlen = 0;
    ireads = 0;
    treads = 0;
    FPRINTF(dbfile,DB_BDATA,0,0)
    if (ALL)
      for (i = 0; i < nreads; i++)
        { rlen = reads[i].rlen;
          if ((reads[i].flags & DB_CCS) == 0)
            css = 0;
          if (rlen >= CUTOFF)
            { if (css == 0 && totlen >= SIZE)
                { FPRINTF(dbfile,DB_BDATA,i,treads)
                  totlen = 0;
                  ireads = 0;
                  nblock += 1;
                }
              ireads += 1;
              treads += 1;
              totlen += rlen;
              css = 1;
            }
        }
    else
      for (i = 0; i < nreads; i++)
        { rlen = reads[i].rlen;
          if (rlen >= CUTOFF && (reads[i].flags & DB_BEST) != 0)
            { ireads += 1;
              treads += 1;
              totlen += rlen;
              if (totlen >= SIZE)
                { FPRINTF(dbfile,DB_BDATA,i+1,treads)
                  totlen = 0;
                  ireads = 0;
                  nblock += 1;
                }
            }
        }

    if (ireads > 0)
      { FPRINTF(dbfile,DB_BDATA,nreads,treads)
        nblock += 1;
      }
    fno = fileno(dbfile);
    FTELLO(upos,dbfile)
    if (ftruncate(fno,upos) < 0)
      SYSTEM_WRITE_ERROR

    FSEEKO(dbfile,dbpos,SEEK_SET)
    FPRINTF(dbfile,DB_NBLOCK,nblock)

    dbs.cutoff = CUTOFF;
    if (ALL)
      dbs.allarr |= DB_ALL;
    else
      dbs.allarr &= ~DB_ALL;
    dbs.treads = treads;
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
