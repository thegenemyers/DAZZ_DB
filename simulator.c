/*******************************************************************************************
 *
 *  Synthetic DNA shotgun dataset simulator
 *     From a supplied reference genome in the form of a Dazzler .dam, sample reads of
 *     mean length -m from a log-normal length distribution with standard deviation -s,
 *     but ignore reads of length less than -x.  Collect enough reads to cover the genome
 *     -c times.   Introduce -e fraction errors into each read where the ratio of insertions,
 *     deletions, and substitutions are set by defined constants INS_RATE and DEL_RATE
 *     within generate.c.  The fraction -f controls the rate at which reads are picked from
 *     the forward and reverse strands which defaults to 50%.  If -C is set then assume the
 *     scaffolds are circular.
 *
 *     The -r parameter seeds the random number generator for the generation of the genome
 *     so that one can reproducbile produce the same underlying genome to sample from.  If
 *     missing, then the job id of the invocation seeds the generator.  The output is sent
 *     to the standard output (i.e. it is a pipe).  The output is in fasta format (i.e. it is
 *     a UNIX pipe).  The output is in Pacbio .fasta format suitable as input to fasta2DB.
 *
 *     The genome is considered a sequence of *scaffolds* (these are reconstituted from the
 *     Dazzler's internal encoding of a .dam), where the gaps are filled with a random
 *     sequence that follows the base distribution of the contigs of the genome.  The program
 *     then samples these filled in scaffolds for reads.  If the -C optioin is set then the
 *     program assumes each scaffold is a circular sequence.
 *
 *     The -M option requests that the scaffold and coordinates from which each read has
 *     been sampled are written to the indicated file, one line per read, ASCII encoded.
 *     This "map" file essentially tells one where every read belongs in an assembly and
 *     is very useful for debugging and testing purposes.  If a read pair is say b,e then
 *     if b < e the read was sampled from [b,e] in the forward direction, and from [e,b]
 *     in the reverse direction otherwise.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *  Mod   :  April 2016 (generates reads w.r.t. a reference genome)
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define PACBIO

#include "DB.h"

static char *Usage[] = { "<genome:dam> [-CU] [-m<int(10000)>]  [-s<int(2000)>] [-e<double(.15)>]",
                         "                   [-c<double(50.)>] [-f<double(.5)>] [-x<int(4000)>]",
                         "                   [-w<int(80)>] [-r<int>] [-M<file>]",
                       };

static int    CIRCULAR;   // -C option
static int    UPPER;      // -U option
static int    RMEAN;      // -m option
static int    RSDEV;      // -s option
static double ERROR;      // -e option
static double COVERAGE;   // -c option
static double FLIP_RATE;  // -f option
static int    RSHORT;     // -x option
static int    WIDTH;      // -w option
static int    HASR;       // -r option is set?
static int    SEED;       // -r option
static FILE  *MAP;        // -M option

#ifdef PACBIO

#define INS_RATE  .73333  // insert rate (for PB data)
#define DEL_RATE  .20000  // deletion rate
#define IDL_RATE  .93333  // insert + delete rate

#elif ILLUMINA

#define INS_RATE  .1  // insert rate (for Illumina data)
#define DEL_RATE  .1  // deletion rate
#define IDL_RATE  .2  // insert + delete rate

#else

#define INS_RATE  .33333  // insert rate (equal weighting)
#define DEL_RATE  .33333  // deletion rate
#define IDL_RATE  .66666  // insert + delete rate

#endif

//  Complement (in the DNA sense) string *s*.

static void complement(int elen, char *s)
{ char *t;
  int   c;

  t = s + (elen-1);
  while (s <= t)
    { c = *s;
      *s = (char) (3-*t);
      *t = (char) (3-c);
      s += 1;
      t -= 1;
    }
}

//  A unit normal distribution random number generator

#define UNORM_LEN 60000
#define UNORM_MAX   6.0

static double unorm_table[UNORM_LEN+1];  // Upper half of cdf of N(0,1)
static double unorm_scale;

static void init_unorm()
{ double del, sum, x;
  int    i;

  unorm_scale = del = UNORM_MAX / UNORM_LEN;

  sum = 0;                            // Integrate pdf, x >= 0 half only.
  for (i = 0; i < UNORM_LEN; i++)
    { x = i * del;
      unorm_table[i] = sum;
      sum += exp(-.5*x*x) * del;
    }
  unorm_table[UNORM_LEN] = sum;

                /* Normalize cdf */
  sum *= 2.;
  for (i = 0; i < UNORM_LEN; i++)
    unorm_table[i] /= sum;
  unorm_table[UNORM_LEN] = 1.;

#ifdef DEBUG
  printf("Truncated tail is < %g\n",
          exp(-.5*UNORM_MAX*UNORM_MAX)/(sum*(1.-exp(-UNORM_MAX))) );
  printf("Diff between last two entries is %g\n",.5-unorm_table[UNORM_LEN-1]);

  printf("\n  CDF:\n");
  for (i = 0; i <= UNORM_LEN; i += 100)
    printf("%6.2f: %10.9f\n",i*del,unorm_table[i]);
#endif
}

static int bin_search(int len, double *tab, double y)
{ int l, m, r;

  // Searches tab[0..len] for min { r : y < tab[r] }.
  //   Assumes y < 1, tab[0] = 0 and tab[len] = 1.
  //   So returned index is in [1,len].

  l = 0;
  r = len;
  while (l < r)
    { m = (l+r) >> 1;
      if (y < tab[m])
        r = m;
      else
        l = m+1;
    }
  return (r);
}

static double sample_unorm(double x)
{ double y;
  int    f;

  if (x >= .5)  // Map [0,1) random var to upper-half of cdf */
    y = x-.5;
  else
    y = .5-x;

  f = bin_search(UNORM_LEN,unorm_table,y);    // Bin. search upper-half cdf
#ifdef DEBUG
  printf("Normal search %g -> %g -> %d",x,y,f);
#endif

  // Linear interpolate between table points

  y = (f - (unorm_table[f]-y) / (unorm_table[f] - unorm_table[f-1]) ) * unorm_scale;

  if (x < .5) y = -y;       // Map upper-half var back to full range
#ifdef DEBUG
  printf(" -> %g\n",y);
#endif

  return (y);
}

//  Open and trim the reference genome *name*.  Determine the number of scaffolds and sizes
//    of each scaffold (in nscaffs and the .coff field of the read records) in the dam.  Then
//    create a sequence for each scaffold (index in the .boff field of the read records), that
//    consists of its contigs with a random sequence filling the gaps (generated according to
//    the bp frequency in db.freq[4]).

DAZZ_DB *load_and_fill(char *name, int *pscaffs)
{ static DAZZ_DB db;
  DAZZ_READ *reads;
  FILE      *bases;
  char      *bases_name;
  char      *seq;
  int        nreads, nscaffs;
  int        i, c;
  int64      ctot;
  int64      o, u;
  double     PRA, PRC, PRG;

  if (Open_DB(name,&db) != 1)
    { fprintf(stderr,"%s: %s is not a Dazzler .dam\n",Prog_Name,name);
      exit (1);
    }
  Trim_DB(&db);

  PRA = db.freq[0];
  PRC = PRA + db.freq[1];
  PRG = PRC + db.freq[2];

  nreads  = db.nreads;
  reads   = db.reads;

  nscaffs = 0;
  for (i = 0; i < nreads; i++)
    if (reads[i].origin == 0)
      nscaffs += 1;

  for (i = 0; i < nscaffs; i++)
    reads[i].coff = 0;

  c = -1;
  for (i = 0; i < nreads; i++)
    { if (reads[i].origin == 0)
        c += 1;
      reads[c].coff = reads[i].fpulse+reads[i].rlen;
    }

  ctot = 0;
  for (i = 0; i < nscaffs; i++)
    ctot += reads[i].coff+1;

  bases_name = Strdup(Catenate(db.path,"","",".bps"),"Allocating base-pair file name");
  bases = Fopen(bases_name,"r");
  if (bases_name == NULL || bases == NULL)
    exit (1);

  seq = (char *) Malloc(ctot+4,"Allocating space for genome");
  if (seq == NULL)
    exit (1);
  *seq++ = 4;

  c = -1;
  o = u = 0;
  for (i = 0; i < nreads; i++)
    { int   len, clen;
      int64 off;

      if (reads[i].origin == 0)
        { if (c >= 0)
            o += reads[c].coff + 1;
          c += 1;
          u = o;
        }
      else
        { int64  p;
          double x;

          p = u + reads[i-1].rlen;
          u = o + reads[i].fpulse;
          while (p < u)
            { x = drand48();
              if (x < PRC)
                if (x < PRA)
                  seq[p++] = 0;
                else
                  seq[p++] = 1;
              else
                if (x < PRG)
                  seq[p++] = 2;
                else
                  seq[p++] = 3;
            }
        }

      len = reads[i].rlen;
      off = reads[i].boff;
      if (ftello(bases) != off)
        FSEEKO(bases,off,SEEK_SET)
      clen = COMPRESSED_LEN(len);
      if (clen > 0)
        FFREAD(seq+u,clen,1,bases)
      Uncompress_Read(len,seq+u);
      if (reads[i].origin == 0)
        reads[c].boff = o;
    }
  reads[nscaffs].boff = ctot;

  db.bases  = (void *) seq;
  db.loaded = 1;

  *pscaffs = nscaffs;
  return (&db);
}

//  Generate reads (a) whose lengths are exponentially distributed with mean *mean* and
//    standard deviation *stdev*, and (b) that are never shorter than *shortest*.  Each
//    read is a randomly sampled interval of one of the filled scaffolds of *source*
//    (each interval is equally likely) that has insertion, deletion, and/or substitution
//    errors introduced into it and which is oriented in either the forward or reverse
//    strand direction with probability FLIP_RATE.  The number of errors introduced is the
//    length of the string times *erate*, and the probability of an insertion, delection,
//    or substitution is controlled by the defined constants INS_RATE and DEL_RATE.
//    If the -C option is set then each scaffold is assumed to be circular and reads can
//    be sampled that span the origin.   Reads are generated until the sum of the lengths of
//    the reads is greater thant coverage times the sum of the lengths of the scaffolds in
//    the reference (i.e. including filled scaffold gaps in the genome size).  The reads are
//    output as fasta entries with the PacBio-specific header format that contains the
//    sampling interval, read length, and a read id.

static void shotgun(DAZZ_DB *source, int nscaffs)
{ DAZZ_READ *reads;
  int64      gleng;
  int        maxlen, nreads, qv;
  int64      totlen, totbp;
  char      *rbuffer, *bases;
  double     nmean, nsdev;
  double    *weights;
  int        scf;

  nsdev = (1.*RSDEV)/RMEAN;
  nsdev = log(1.+nsdev*nsdev);
  nmean = log(1.*RMEAN) - .5*nsdev;
  nsdev = sqrt(nsdev);

  bases = source->bases;
  reads = source->reads;
  gleng = reads[nscaffs].boff - nscaffs;
  if (gleng <= RSHORT)
    { fprintf(stderr,"Genome length is less than shortest read length !\n");
      exit (1);
    }

  init_unorm();

  weights = (double *) Malloc(sizeof(double)*(nscaffs+1),"Allocating contig weights");
  if (weights == NULL)
    exit (1);

  { double r;

    r = 0.;
    for (scf = 0; scf < nscaffs; scf++)
      { weights[scf] = r/gleng;
        r += reads[scf].coff;
      }
    weights[nscaffs] = 1.;
  }

  qv = (int) (1000 * (1.-ERROR));

  rbuffer = NULL;
  maxlen  = 0;
  totlen  = 0;
  totbp   = COVERAGE*gleng;
  nreads  = 0;
  while (totlen < totbp)
    { int    len, sdl, ins, del, elen, slen, rbeg, rend;
      int    j;
      double uni;
      char  *s, *t;

      scf = bin_search(nscaffs,weights,drand48()) - 1;   //  Pick a scaffold with probabilitye
                                                         //    proportional to its length

      uni = drand48();
      len = (int) exp(nmean + nsdev*sample_unorm(uni));    //  Pick a read length
      if (len <= RSHORT)
        continue;

      // New sampler:

      slen = reads[scf].coff;
      rbeg = (int) (drand48()*slen);          //  Pick a spot for read start
      if (CIRCULAR)
        rend = (rbeg + len) % slen;           //  Wrap if circular
      else
        { if (drand48() < .5)                 //  Pick direction and trim if necessary
            { rend = rbeg + len;              //    if not circular
              if (rend > slen)
                { rend = slen;
                  len  = rend - rbeg;
                }
            }
          else
            { rend = rbeg;
              rbeg = rbeg - len;
              if (rbeg < 0)
                { rbeg = 0;
                  len  = rend;
                }
            }
          if (len <= RSHORT)
            continue;
        }

      // Old sampler:
      //
      // rbeg = (int) (drand48()*((reads[scf].coff-len)+.9999999));
      // rend = rbeg + len;

      sdl = (int) (len*ERROR);      //  Determine number of inserts *ins*, deletions *del,
      ins = del = 0;                //    and substitions+deletions *sdl*.
      for (j = 0; j < sdl; j++)
        { double x = drand48();
          if (x < INS_RATE)
            ins += 1;
          else if (x < IDL_RATE)
            del += 1;
        }
      sdl -= ins;
      elen = len + (ins-del);

      if (elen > maxlen)
        { maxlen  = ((int) (1.2*elen)) + 1000;
          rbuffer = (char *) Realloc(rbuffer,maxlen+3,"Allocating read buffer");
          if (rbuffer == NULL)
            exit (1);
        }

      t = rbuffer;
      s = bases + (reads[scf].boff + rbeg);

      //   Generate the string with errors.  NB that inserts occur randomly between source
      //     characters, while deletions and substitutions occur on source characters.

      while ((len+1) * drand48() < ins)
        { *t++ = (char) (4.*drand48());
          ins -= 1;
        }
      for ( ; len > 0; len--)
        { if (len * drand48() >= sdl)
            *t++ = *s;
          else if (sdl * drand48() >= del)
            { double x = 3.*drand48();
              if (x >= *s)
                x += 1.;
              *t++ = (char) x;
              sdl -= 1;
            }
          else
            { del -= 1;
              sdl -= 1;
            }
          s += 1;
          if (*s == 4)
            s = bases + reads[scf].boff;
          while (len * drand48() < ins)
            { *t++ = (char) (4.*drand48());
              ins -= 1;
            }
        }
      *t = 4;

      if (drand48() >= FLIP_RATE)    //  Complement the string with probability FLIP_RATE.
        { complement(elen,rbuffer);
          j = rend;
          rend = rbeg;
          rbeg = j;
        }

      PRINTF(">Sim/%d/%d_%d RQ=0.%d\n",nreads+1,0,elen,qv)
      if (UPPER)
        Upper_Read(rbuffer);
      else
        Lower_Read(rbuffer);
      for (j = 0; j+WIDTH < elen; j += WIDTH)
        PRINTF("%.*s\n",WIDTH,rbuffer+j)
      if (j < elen)
        PRINTF("%s\n",rbuffer+j)

       if (MAP != NULL)
         FPRINTF(MAP," %6d %9d %9d\n",scf,rbeg,rend)

       totlen += elen;
       nreads += 1;
    }
}

int main(int argc, char *argv[])
{ DAZZ_DB *source;
  int      nscaffs;

  //  Process command line

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("simulator");

    RMEAN     = 10000;
    RSDEV     = 2000;
    ERROR     = .15;
    COVERAGE  = 50.;
    FLIP_RATE = .5;
    RSHORT    = 4000;
    HASR      = 0;
    MAP       = NULL;
    WIDTH     = 80;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("CU");
            break;
          case 'c':
            ARG_REAL(COVERAGE)
            if (COVERAGE < 0.)
              { fprintf(stderr,"%s: Coverage must be non-negative (%g)\n",Prog_Name,COVERAGE);
                exit (1);
              }
            break;
          case 'e':
            ARG_REAL(ERROR)
            if (ERROR < 0. || ERROR > .5)
              { fprintf(stderr,"%s: Error rate must be in [0,.5] (%g)\n",Prog_Name,ERROR);
                exit (1);
              }
            break;
          case 'f':
            ARG_REAL(FLIP_RATE)
            if (FLIP_RATE < 0. || FLIP_RATE > 1.)
              { fprintf(stderr,"%s: Error rate must be in [0,1] (%g)\n",Prog_Name,FLIP_RATE);
                exit (1);
              }
            break;
          case 'm':
            ARG_POSITIVE(RMEAN,"Mean read length")
            break;
          case 'r':
            SEED = strtol(argv[i]+2,&eptr,10);
            HASR = 1;
            if (*eptr != '\0' || argv[i][2] == '\0')
              { fprintf(stderr,"%s: -r argument is not an integer\n",Prog_Name);
                exit (1);
              }
            break;
          case 's':
            ARG_NON_NEGATIVE(RSDEV,"Read length standard deviation")
            break;
          case 'x':
            ARG_NON_NEGATIVE(RSHORT,"Read length minimum")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
          case 'M':
            MAP = Fopen(argv[i]+2,"w");
            if (MAP == NULL)
              exit (1);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    CIRCULAR = flags['C'];
    UPPER    = flags['U'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -m: average read length (log normal distribution).\n");
        fprintf(stderr,"      -s: standard deviation of read lengths (log normal)\n");
        fprintf(stderr,"      -x: ignore reads below this length\n");
        fprintf(stderr,"      -f: forward/reverse strand sampling fraction\n");
        fprintf(stderr,"      -e: error rate\n");
        fprintf(stderr,"      -c: coverage of genome\n");
        fprintf(stderr,"      -C: assume genome is circular (default is linear)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -r: Random number generator seed (default is process id).\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        fprintf(stderr,"      -M: create a map file that indicates where every read was sampled\n");
        exit (1);
      }
  }

  if (HASR)
    srand48(SEED);
  else
    srand48(getpid());

  //  Read and generate

  source = load_and_fill(argv[1],&nscaffs);

  shotgun(source,nscaffs);

  if (MAP != NULL)
    fclose(MAP);

  exit (0);
}
