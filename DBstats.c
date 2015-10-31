/************************************************************************************\
*                                                                                    *
* Copyright (c) 2014, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

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

static char *Usage = " [-nu] [-b<int(1000)>] [-m<track>]+ <name:db|dam>";

int main(int argc, char *argv[])
{ HITS_DB _db, *db = &_db;
  int     dam;

  int64   ototal;
  int     oreads;
  int     nbin, *hist;
  int64  *bsum;

  int     NONE;
  int     TRIM;
  int     BIN;

  int     MMAX, MTOP;
  char  **MASK;

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("DBstats")

    BIN    = 1000;
    MTOP  = 0;
    MMAX  = 10;
    MASK  = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    if (MASK == NULL)
      exit (1);

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("nu")
            break;
          case 'b':
            ARG_POSITIVE(BIN,"Bin size")
            break;
          case 'm':
            if (MTOP >= MMAX)
              { MMAX = 1.2*MTOP + 10;
                MASK = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
                if (MASK == NULL)
                  exit (1);
              }
            MASK[MTOP++] = argv[i]+2;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    NONE = flags['n'];
    TRIM = 1-flags['u'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { int i, status, kind;

    //  Open .db or .dam

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    dam = status;

    //  Check tracks and load tracks for untrimmed DB

    for (i = 0; i < MTOP; i++)
      { status = Check_Track(db,MASK[i],&kind);
        if (status == -2)
          fprintf(stderr,"%s: Warning: -m%s option given but no track found.\n",Prog_Name,MASK[i]);
        else if (status == -1)
          fprintf(stderr,"%s: Warning: %s track not sync'd with db.\n",Prog_Name,MASK[i]);
        else if (kind != MASK_TRACK)
          fprintf(stderr,"%s: Warning: %s track is not a mask track.\n",Prog_Name,MASK[i]);
        else if (status == 0)
          Load_Track(db,MASK[i]);
        else if (status == 1 && !TRIM)
          fprintf(stderr,"%s: Warning: %s track is for a trimmed db but -u is set.\n",
                         Prog_Name,MASK[i]);
      }

    oreads = db->nreads;
    ototal = db->totlen;

    if (TRIM)
      { Trim_DB(db);

        //  Load tracks for trimmed DB

        for (i = 0; i < MTOP; i++)
          { status = Check_Track(db,MASK[i],&kind);
            if (status < 0)
              continue;
            else if (status == 1)
              Load_Track(db,MASK[i]);
          }
      }
  }

  { int        i;
    int64      totlen;
    int        nreads, maxlen;
    int64      ave, dev;
    HITS_READ *reads;

    nreads = db->nreads;
    totlen = db->totlen;
    maxlen = db->maxlen;
    reads  = db->reads;

    nbin  = (maxlen-1)/BIN + 1;
    hist  = (int *) Malloc(sizeof(int)*nbin,"Allocating histograms");
    bsum  = (int64 *) Malloc(sizeof(int64)*nbin,"Allocating histograms");
    if (hist == NULL || bsum == NULL)
      exit (1);

    for (i = 0; i < nbin; i++)
      { hist[i] = 0;
        bsum[i] = 0;
      }
 
    for (i = 0; i < nreads; i++)
      { int rlen = reads[i].rlen;
        hist[rlen/BIN] += 1;
        bsum[rlen/BIN] += rlen;
      }

    nbin = (maxlen-1)/BIN + 1;
    ave  = totlen/nreads;
    dev  = 0;
    for (i = 0; i < nreads; i++)
      { int rlen = reads[i].rlen;
        dev += (rlen-ave)*(rlen-ave);
      }
    dev = (int64) sqrt((1.*dev)/nreads);

    if (dam)
      printf("\nStatistics for all contigs");
    else if (db->all || !TRIM)
      printf("\nStatistics for all wells");
    else
      printf("\nStatistics for all reads");
    if (TRIM && db->cutoff > 0)
      { printf(" of length ");
        Print_Number(db->cutoff,0,stdout);
        printf(" bases or more\n\n");
      }
    else if (dam)
      printf(" in the map index\n\n");
    else
      printf(" in the data set\n\n");

    Print_Number((int64) nreads,15,stdout);
    if (dam)
      printf(" contigs");
    else
      printf(" reads  ");
    if (TRIM)
      { printf("      out of ");
        Print_Number((int64 ) oreads,15,stdout);
        printf("  (%5.1f%%)",(100.*nreads)/oreads);
      }
    printf("\n");

    Print_Number(totlen,15,stdout);
    printf(" base pairs");
    if (TRIM)
      { printf("   out of ");
        Print_Number(ototal,15,stdout);
        printf("  (%5.1f%%)",(100.*totlen)/ototal);
      }
    printf("\n\n");

    Print_Number(ave,15,stdout);
    if (dam)
      printf(" average contig length\n");
    else
      { printf(" average read length\n");
        Print_Number(dev,15,stdout);
        printf(" standard deviation\n");
      }

    printf("\n  Base composition: %.3f(A) %.3f(C) %.3f(G) %.3f(T)\n",
           db->freq[0],db->freq[1],db->freq[2],db->freq[3]);

    if (!NONE)
      { int64 btot;
        int   cum, skip;

        printf("\n  Distribution of Read Lengths (Bin size = ");
        Print_Number((int64) BIN,0,stdout);
        printf(")\n\n        Bin:      Count  %% Reads  %% Bases     Average\n");
        if (dam)
          skip = 0;
        else
          skip = -1;
        cum  = 0;
        btot = 0;
        for (i = nbin-1; i >= 0; i--)
          { cum  += hist[i];
            btot += bsum[i];
            if (hist[i] != skip)
              { Print_Number((int64) (i*BIN),11,stdout);
                printf(":");
                Print_Number((int64) hist[i],11,stdout);
                printf("    %5.1f    %5.1f   %9lld\n",(100.*cum)/nreads,
                                                      (100.*btot)/totlen,btot/cum);
              }
            if (cum == nreads) break;
          }
      }
  }

  { int64      totlen;
    int        numint, maxlen;
    int64      ave, dev;
    HITS_TRACK *track;

    for (track = db->tracks; track != NULL; track = track->next)
      { char  *data = track->data;
        int64 *anno = (int64 *) track->anno;
        int    k, rlen;
        int   *idata, *edata;

        totlen = 0;
        numint = 0;
        maxlen = 0;
        for (k = 0; k < db->nreads; k++)
          { edata = (int *) (data + anno[k+1]);
            for (idata = (int *) (data + anno[k]); idata < edata; idata += 2)
              { rlen = idata[1] - *idata;
                numint += 1;
                totlen += rlen;
                if (rlen > maxlen)
                  maxlen = rlen;
              }
          }

        nbin = (maxlen-1)/BIN + 1;

        for (k = 0; k < nbin; k++)
          { hist[k] = 0;
            bsum[k] = 0;
          }

        ave  = totlen/numint;
        dev  = 0;
        for (k = 0; k < db->nreads; k++)
          { edata = (int *) (data + anno[k+1]);
            for (idata = (int *) (data + anno[k]); idata < edata; idata += 2)
              { rlen = idata[1] - *idata;
                dev += (rlen-ave)*(rlen-ave);
                hist[rlen/BIN] += 1;
                bsum[rlen/BIN] += rlen;
              }
          }
        dev = (int64) sqrt((1.*dev)/numint);

        printf("\n\nStatistics for %s-track\n",track->name);

        printf("\n  There are ");
        Print_Number(numint,0,stdout);
        printf(" intervals totaling ");
        Print_Number(totlen,0,stdout);
        printf(" bases (%.1f%% of all data)\n",(100.*totlen)/db->totlen);

        { int64 btot;
          int   cum;

          printf("\n  Distribution of %s intervals (Bin size = ",track->name);
          Print_Number((int64) BIN,0,stdout);
          printf(")\n\n        Bin:      Count  %% Intervals  %% Bases     Average\n");
          cum  = 0;
          btot = 0;
          for (k = nbin-1; k >= 0; k--)
            { cum  += hist[k];
              btot += bsum[k];
              if (hist[k] > 0)
                { Print_Number((int64) (k*BIN),11,stdout);
                  printf(":");
                  Print_Number((int64) hist[k],11,stdout);
                  printf("        %5.1f    %5.1f   %9lld\n",(100.*cum)/numint,
                                                        (100.*btot)/totlen,btot/cum);
                  if (cum == numint) break;
                }
            }
          printf("\n");
        }
      }
  }

  free(hist);
  free(bsum);
  Close_DB(db);

  exit (0);
}
