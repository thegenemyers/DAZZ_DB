/*******************************************************************************************
 *
 *  My implementation of the SDUST algorithm (Morgulis et al., JCB 13, 5 (2006), 1028-1040)
 *
 *  Author:  Gene Myers
 *  Date  :  September 2013
 *  Mod   :  Is now incremental
 *  Date  :  April 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

#undef DEBUG

#ifdef DEBUG

static int Caps[4] = { 'A', 'C', 'G', 'T' };
static int Lowr[4] = { 'a', 'c', 'g', 't' };

#endif

static char *Usage = "[-b] [-w<int(64)>] [-t<double(2.)>] [-m<int(10)>] <path:db|dam>";

typedef struct _cand
  { struct _cand *next;
    struct _cand *prev;
    int           beg;
    int           end;
    double        score;
  } Candidate;

int main(int argc, char *argv[])
{ DAZZ_DB   _db, *db = &_db;
  FILE      *afile, *dfile;
  int64      indx;
  int        nreads;
  int       *mask;
  Candidate *cptr;

  int        WINDOW;
  double     THRESH;
  int        MINLEN;
  int        BIASED;

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("DBdust")

    WINDOW = 64;
    THRESH = 2.;
    MINLEN = 9;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("b")
            break;
          case 'w':
            ARG_POSITIVE(WINDOW,"Window size")
            break;
          case 't':
            ARG_REAL(THRESH)
            if (THRESH <= 0.)
              { fprintf(stderr,"%s: Threshold must be positive (%g)\n",Prog_Name,THRESH);
                exit (1);
              }
            break;
          case 'm':
            ARG_NON_NEGATIVE(MINLEN,"Minimum hit")
            MINLEN -= 1;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    BIASED = flags['b'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -w: DUST algorithm window size.\n");
        fprintf(stderr,"      -t: DUST algorithm threshold.\n");
        fprintf(stderr,"      -m: Record only low-complexity intervals >= this size.\n");
        fprintf(stderr,"      -b: Take into account base composition bias.\n");
        exit (1);
      }
  }

  //  Open .db or .dam

  { int status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
  }

  mask = (int *) Malloc((db->maxlen+1)*sizeof(int),"Allocating mask vector");
  cptr = (Candidate *) Malloc((WINDOW+1)*sizeof(Candidate),"Allocating candidate vector");
  if (mask == NULL || cptr == NULL)
    exit (1);

  { char *pwd, *root, *fname;
    int   size;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".db");
    size  = 0;

    fname = Catenate(pwd,PATHSEP,root,".dust.anno");
    if ((afile = fopen(fname,"r+")) == NULL || db->part > 0)
      { if (afile != NULL)
          fclose(afile);
        afile = Fopen(fname,"w");
        dfile = Fopen(Catenate(pwd,PATHSEP,root,".dust.data"),"w");
        if (dfile == NULL || afile == NULL)
          exit (1);
        FFWRITE(&(db->nreads),sizeof(int),1,afile)
        FFWRITE(&size,sizeof(int),1,afile)
        nreads = 0;
        indx = 0;
        FFWRITE(&indx,sizeof(int64),1,afile)
      }
    else
      { dfile = Fopen(Catenate(pwd,PATHSEP,root,".dust.data"),"r+");
        if (dfile == NULL)
          exit (1);
        if (fread(&nreads,sizeof(int),1,afile) != 1)
          SYSTEM_READ_ERROR
        if (nreads >= db->nreads)
          { fclose(afile);
            fclose(dfile);
            exit(0);
          }
        FSEEKO(afile,0,SEEK_SET)
        FFWRITE(&(db->nreads),sizeof(int),1,afile)
        FFWRITE(&size,sizeof(int),1,afile)
        FSEEKO(afile,0,SEEK_END)
        FSEEKO(dfile,0,SEEK_END)
        FTELLO(indx,dfile)
      }

    free(pwd);
    free(root);
  }

  { int       *mask1;
    char      *read, *lag2;
    int        wcount[64], lcount[64];
    Candidate *aptr;
    double     skew[64], thresh2r;
    int        thresh2i;
    int        i;

    read = New_Read_Buffer(db);
    lag2 = read-2;

    mask1 = mask+1;
    *mask = -2;

    aptr  = cptr+1;
    for (i = 1; i < WINDOW; i++)
      cptr[i].next = aptr+i;
    cptr[WINDOW].next = NULL;

    cptr->next = cptr->prev = cptr;
    cptr->beg  = -2;

    thresh2r = 2.*THRESH;
    thresh2i = (int) ceil(thresh2r);

    if (BIASED)
      { int a, b, c, p;

        p = 0;
        for (a = 0; a < 4; a++)
         for (b = 0; b < 4; b++)
          for (c = 0; c < 4; c++)
            skew[p++] = .015625 / (db->freq[a]*db->freq[b]*db->freq[c]);
      }

    for (i = nreads; i < db->nreads; i++)
      { Candidate *lptr, *jptr;
        int       *mtop;
        double     mscore;
        int        len;
        int        wb, lb;
        int        j, c, d;

        len = db->reads[i].rlen;	  //  Fetch read
        Load_Read(db,i,read,0);

        c = (read[0] << 2) | read[1];     //   Convert to triple codes
        for (j = 2; j < len; j++)
          { c = ((c << 2) & 0x3f) | read[j];
            lag2[j] = (char) c;
          }
        len -= 2;

        for (j = 0; j < 64; j++)		//   Setup counter arrays
          wcount[j] = lcount[j] = 0;

        mtop = mask;                      //   The dust algorithm
        lb   = wb   = -1;

        if (BIASED)

          { double lsqr, wsqr, trun;      //   Modification for high-compositional bias

            wsqr = lsqr = 0.;
            for (j = 0; j < len; j++)
              { c = read[j];

#define ADDR(e,cnt,sqr)	 sqr += (cnt[e]++) * skew[e];

#define DELR(e,cnt,sqr)	 sqr -= (--cnt[e]) * skew[e];

#define WADDR(e) ADDR(e,wcount,wsqr)
#define WDELR(e) DELR(e,wcount,wsqr)
#define LADDR(e) ADDR(e,lcount,lsqr)
#define LDELR(e) DELR(e,lcount,lsqr)

                if (j > WINDOW-3)
                  { d = read[++wb];
                    WDELR(d)
                  }
                WADDR(c)

                if (lb < wb)
                  { d = read[++lb];
                    LDELR(d)
                  }
                trun  = (lcount[c]++) * skew[c];
                lsqr += trun;
                if (trun >= thresh2r)
                  { while (lb < j)
                      { d = read[++lb];
                        LDELR(d)
                        if (d == c) break;
                      }
                  }

                jptr = cptr->prev;
                if (jptr != cptr && jptr->beg <= wb)
                  { c = jptr->end + 2;
                    if (*mtop+1 >= jptr->beg)
                      { if (*mtop < c)
                          *mtop = c;
                      }
                    else
                      { *++mtop = jptr->beg;
                        *++mtop = c;
                      }
                    lptr = jptr->prev;
                    cptr->prev = lptr;
                    lptr->next = cptr;
                    jptr->next = aptr;
                    aptr = jptr;
                  }

                if (wsqr <= lsqr*THRESH) continue;

                jptr   = cptr->next;
                lptr   = cptr;
                mscore = 0.;
                for (c = lb; c > wb; c--)
                  { d = read[c];
                    LADDR(d)
                    if (lsqr >= THRESH * (j-c))
                      { for ( ; jptr->beg >= c; jptr = (lptr = jptr)->next)
                          if (jptr->score > mscore)
                            mscore = jptr->score;
                        if (lsqr >= mscore * (j-c))
                          { mscore = lsqr / (j-c);
                            if (lptr->beg == c)
                              { lptr->end   = j;
                                lptr->score = mscore;
                              }
                            else
                              { aptr->beg   = c;
                                aptr->end   = j;
                                aptr->score = mscore;
                                aptr->prev  = lptr;
                                lptr = lptr->next = aptr;
                                aptr = aptr->next;
                                jptr->prev = lptr;
                                lptr->next = jptr;
                              }
                          }
                      }
                  }

                for (c++; c <= lb; c++)
                  { d = read[c];
                    LDELR(d)
                  }
              }
          }

        else

          { int lsqr, wsqr, trun;                 //  Algorithm for GC-balanced sequences

            wsqr = lsqr = 0;
            for (j = 0; j < len; j++)
              { c = read[j];

#define ADDI(e,cnt,sqr)	 sqr += (cnt[e]++);

#define DELI(e,cnt,sqr)	 sqr -= (--cnt[e]);

#define WADDI(e) ADDI(e,wcount,wsqr)
#define WDELI(e) DELI(e,wcount,wsqr)
#define LADDI(e) ADDI(e,lcount,lsqr)
#define LDELI(e) DELI(e,lcount,lsqr)

                if (j > WINDOW-3)
                  { d = read[++wb];
                    WDELI(d)
                  }
                WADDI(c)

                if (lb < wb)
                  { d = read[++lb];
                    LDELI(d)
                  }
                trun  = lcount[c]++;
                lsqr += trun;
                if (trun >= thresh2i)
                  { while (lb < j)
                      { d = read[++lb];
                        LDELI(d)
                        if (d == c) break;
                      }
                  }

                jptr = cptr->prev;
                if (jptr != cptr && jptr->beg <= wb)
                  { c = jptr->end + 2;
                    if (*mtop+1 >= jptr->beg)
                      { if (*mtop < c)
                          *mtop = c;
                      }
                    else
                      { *++mtop = jptr->beg;
                        *++mtop = c;
                      }
                    lptr = jptr->prev;
                    cptr->prev = lptr;
                    lptr->next = cptr;
                    jptr->next = aptr;
                    aptr = jptr;
                  }

                if (wsqr <= lsqr*THRESH) continue;

                jptr   = cptr->next;
                lptr   = cptr;
                mscore = 0.;
                for (c = lb; c > wb; c--)
                  { d = read[c];
                    LADDI(d)
                    if (lsqr >= THRESH * (j-c))
                      { for ( ; jptr->beg >= c; jptr = (lptr = jptr)->next)
                          if (jptr->score > mscore)
                            mscore = jptr->score;
                        if (lsqr >= mscore * (j-c))
                          { mscore = (1. * lsqr) / (j-c);
                            if (lptr->beg == c)
                              { lptr->end   = j;
                                lptr->score = mscore;
                              }
                            else
                              { aptr->beg   = c;
                                aptr->end   = j;
                                aptr->score = mscore;
                                aptr->prev  = lptr;
                                lptr = lptr->next = aptr;
                                aptr = aptr->next;
                                jptr->prev = lptr;
                                lptr->next = jptr;
                              }
                          }
                      }
                  }

                for (c++; c <= lb; c++)
                  { d = read[c];
                    LDELI(d)
                  }
              }
          }

        while ((jptr = cptr->prev) != cptr)
          { c = jptr->end + 2;
            if (*mtop+1 >= jptr->beg)
              { if (*mtop < c)
                  *mtop = c;
              }
            else
              { *++mtop = jptr->beg;
                *++mtop = c;
              }
            cptr->prev = jptr->prev;
            jptr->prev->next = cptr;
            jptr->next = aptr;
            aptr = jptr;
          }

        { int *jtop, ntop;

          ntop = 0;
          for (jtop = mask1; jtop < mtop; jtop += 2)
            if (jtop[1] - jtop[0] >= MINLEN)
              { mask[++ntop] = jtop[0];
                mask[++ntop] = jtop[1]+1;
              }
          mtop  = mask + ntop;
          indx += ntop*sizeof(int);
          FFWRITE(&indx,sizeof(int64),1,afile)
          FFWRITE(mask1,sizeof(int),ntop,dfile)
        }

#ifdef DEBUG

        { int *jtop;

          printf("\nREAD %d\n",i);
          for (jtop = mask1; jtop < mtop; jtop += 2)
            printf(" [%5d,%5d]\n",jtop[0],jtop[1]);

          Load_Read(db,i,read,0);

          jtop = mask1;
          for (c = 0; c < len; c++)
            { while (jtop < mtop && c > jtop[1])
                jtop += 2;
              if (jtop < mtop && c >= *jtop)
                printf("%c",Caps[(int) read[c]]);
              else
                printf("%c",Lowr[(int) read[c]]);
              if ((c%80) == 79)
                printf("\n");
            }
          printf("\n");
        }

#endif
      }
  }

  FCLOSE(afile)
  FCLOSE(dfile)

  Close_DB(db);

  exit (0);
}
