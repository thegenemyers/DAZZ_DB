/*******************************************************************************************
 *
 *  Display statistics about the contents of a .db and a histogram of its read lengths.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *  Mod   :  April 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "DB.h"

static char *Usage = " [-a] [-x<int>] [-b<int(1000)>] <name:db>";

int main(int argc, char *argv[])
{ HITS_DB     db;
  HITS_READ  *reads;

  int        ALL;
  int        CUTOFF;
  int        BIN;

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("DBstats")

    CUTOFF = 0;
    BIN    = 1000;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("a")
            break;
          case 'x':
            ARG_NON_NEGATIVE(CUTOFF,"Min read length cutoff")
            break;
          case 'b':
            ARG_POSITIVE(BIN,"Bin size")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    ALL = flags['a'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  if (Open_DB(argv[1],&db))
    exit (1);
  reads = db.reads;

  { int         nbin, *hist;
    int64       totlen, *bsum;
    int         nreads;
    int         i;

    nbin = (db.maxlen-1)/BIN + 1;
    hist = (int *) Malloc(sizeof(int)*nbin,"Allocating histograms");
    bsum = (int64 *) Malloc(sizeof(int64)*nbin,"Allocating histograms");
    if (hist == NULL || bsum == NULL)
      exit (1);

    for (i = 0; i < nbin; i++)
      { hist[i] = 0;
        bsum[i] = 0;
      }

    { int64 ave, dev;
 
      if ( ! ALL)
        { for (i = 0; i < db.nreads;  i++)
            if ((reads[i].flags & DB_BEST) == 0)
              { reads[i].beg = 1;
                reads[i].end = 0;
              }
        }

      totlen = 0;
      nreads = 0;
      for (i = 0; i < db.nreads; i++)
        { int rlen = ((int) reads[i].end) - reads[i].beg;
          if (rlen >= CUTOFF)
            { totlen += rlen;
              nreads += 1;
              hist[rlen/BIN] += 1;
              bsum[rlen/BIN] += rlen;
            }
        }

      ave = totlen/nreads;
      dev = 0;
      for (i = 0; i < db.nreads; i++)
        { int rlen = ((int) reads[i].end) - reads[i].beg;
          if (rlen >= CUTOFF)
            dev += (rlen-ave)*(rlen-ave);
        }
      dev = sqrt(dev/nreads);

      if (CUTOFF == 0 && ALL)
        { printf("\nStatistics over all reads in the data set\n\n");
          Print_Number((int64) nreads,15);
          printf(" reads\n");
          Print_Number(totlen,15);
          printf(" base pairs\n");
        }
      else
        { if (ALL)
            printf("\nStatistics for all reads");
          else
            printf("\nStatistics for all wells");
          if (CUTOFF > 0)
            { printf(" of length ");
              Print_Number(CUTOFF,0);
              printf(" bases or more\n\n");
            }
          else
            printf(" in the data set\n\n");
          Print_Number((int64) nreads,15);
          printf(" reads       out of ");
          Print_Number((int64 ) db.nreads,15);
          printf("   %5.1f%% loss\n",100.-(100.*nreads)/db.nreads);
          Print_Number(totlen,15);
          printf(" base pairs  out of ");
          Print_Number(db.totlen,15);
          printf("   %5.1f%% loss\n",100.-(100.*totlen)/db.totlen);
        }
      printf("\nBase composition: %.3f(A) %.3f(C) %.3f(G) %.3f(T)\n\n",
             db.freq[0],db.freq[1],db.freq[2],db.freq[3]);
      Print_Number(ave,15);
      printf(" average read length\n");
      Print_Number(dev,15);
      printf(" standard deviation\n");
    }

    { int64 btot;
      int   cum;

      printf("\nDistribution of Read Lengths (Bin size = ");
      Print_Number((int64) BIN,0);
      printf(")\n\n    Bin:      Count  %% Reads  %% Bases   Average\n");
      cum  = 0;
      btot = 0;
      for (i = nbin-1; i >= 0; i--)
        { cum  += hist[i];
          btot += bsum[i];
          Print_Number((int64) (i*BIN),7);
          printf(":");
          Print_Number((int64) hist[i],11);
          printf("    %5.1f    %5.1f    %5lld\n",(100.*cum)/nreads,(100.*btot)/totlen,btot/cum);
          if (cum == nreads) break;
        }
      printf("\n");
    }
  }

  { HITS_TRACK *dust;

    dust = Load_Track(&db,"dust");
    if (dust != NULL)
      { void *data = dust->data;
        int  *anno = (int *) dust->anno;
        int   i, rlen;
        int  *idata, *edata;
        int64 numint, totlen;

        totlen = 0;
        numint = 0;
        for (i = 0; i < db.nreads; i++)
          { rlen = reads[i].end - reads[i].beg;
            if (rlen >= CUTOFF)
              { edata = (int *) (data + anno[i+1]);
                for (idata = (int *) (data + anno[i]); idata < edata; idata += 2)
                  { numint += 1;
                    totlen += (idata[1] - *idata) + 1;
                  }
              }
          }
        Print_Number(numint,0);
        printf(" low-complexity intervals totaling ");
        Print_Number(totlen,0);
        printf(" bases\n\n");
      }
    Close_Track(&db,"dust");
  }

  Close_DB(&db);

  exit (0);
}
