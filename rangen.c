/*******************************************************************************************
 *
 *  Synthetic DNA shotgun sequence generator
 *     Generate a fake genome of size genlen*1Mb long, that has an AT-bias of -b.
 *     The -r parameter seeds the random number generator for the generation of the genome
 *     so that one can reproducbile produce the same underlying genome to sample from.  If
 *     missing, then the job id of the invocation seeds the generator.  The sequence is
 *     sent to the standard output in .fasta format.
 *
 *  Author:  Gene Myers
 *  Date  :  April 2016
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

static char *Usage = "<genlen:double> [-U] [-b<double(.5)>] [-w<int(80)>] [-r<int>]";

static int    GENOME;     // -g option * 1Mbp
static double BIAS;       // -b option
static int    HASR = 0;   // -r option is set?
static int    SEED;       // -r option
static int    WIDTH;      // -w option
static int    UPPER;      // -U option

static char *Prog_Name;

//  Generate a random DNA sequence of length *len* with an AT-bias of BIAS.
//    Uppercase letters if UPPER is set, lowercase otherwise.

static char *random_genome(int len)
{ static char  *seq = NULL;
  static double x, PRA, PRC, PRG;
  int    i;

  if (seq == NULL)
    { PRA = BIAS/2.;
      PRC = (1.-BIAS)/2. + PRA;
      PRG = (1.-BIAS)/2. + PRC;
      if ((seq = (char *) malloc(WIDTH+1)) == NULL)
        { fprintf(stderr,"%s: Allocating genome sequence\n",Prog_Name);
          exit (1);
        }
    }

  if (UPPER)
    for (i = 0; i < len; i++)
      { x = drand48();
        if (x < PRC)
          if (x < PRA)
            seq[i] = 'A';
          else
            seq[i] = 'C';
        else
          if (x < PRG)
            seq[i] = 'G';
          else
            seq[i] = 'T';
      }
  else
    for (i = 0; i < len; i++)
      { x = drand48();
        if (x < PRC)
          if (x < PRA)
            seq[i] = 'a';
          else
            seq[i] = 'c';
        else
          if (x < PRG)
            seq[i] = 'g';
          else
            seq[i] = 't';
      }
  seq[len] = '\0';

  return (seq);
}

int main(int argc, char *argv[])
{ int    i, j;
  char  *eptr;
  double glen;
 
  //  Process command arguments
  //
  //  Usage: <GenomeLen:double> [-b<double(.5)>] [-r<int>]

  Prog_Name = strdup("rangen");

  WIDTH = 80;
  BIAS  = .5;
  HASR  = 0;
  UPPER = 0;

  j = 1;
  for (i = 1; i < argc; i++)
    if (argv[i][0] == '-')
      switch (argv[i][1])
      { default:
          fprintf(stderr,"%s: %s is an illegal option\n",Prog_Name,argv[i]);
          exit (1);
        case 'U':
          if (argv[i][2] != '\0')
            { fprintf(stderr,"%s: %s is an illegal option\n",Prog_Name,argv[i]);
              exit (1);
            }
          UPPER = 1;
          break;
        case 'b':
          BIAS = strtod(argv[i]+2,&eptr);
          if (*eptr != '\0' || argv[i][2] == '\0')
            { fprintf(stderr,"%s: -%c '%s' argument is not a real number\n",
                             Prog_Name,argv[i][1],argv[i]+2);
              exit (1);
            }
          if (BIAS < 0. || BIAS > 1.)
            { fprintf(stderr,"%s: AT-bias must be in [0,1] (%g)\n",Prog_Name,BIAS);
              exit (1);
            }
          break;
        case 'r':
          SEED = strtol(argv[i]+2,&eptr,10);
          HASR = 1;
          if (*eptr != '\0' || argv[i][2] == '\0')
            { fprintf(stderr,"%s: -r argument is not an integer\n",Prog_Name);
              exit (1);
            }
          break;
        case 'w':
          WIDTH = strtol(argv[i]+2,&eptr,10);
          if (*eptr != '\0' || argv[i][2] == '\0')
            { fprintf(stderr,"%s: -w '%s' argument is not an integer\n",Prog_Name,argv[i]+2);
              exit (1);
            }
          if (WIDTH < 0)
            { fprintf(stderr,"%s: Line width must be non-negative (%d)\n",Prog_Name,WIDTH);
              exit (1);
            }
          break;
      }
    else
      argv[j++] = argv[i];
  argc = j;

  if (argc != 2)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      fprintf(stderr,"\n");
      fprintf(stderr,"      -r: Random number generator seed (default is process id).\n");
      fprintf(stderr,"      -b: AT vs GC ratio or bias.\n");
      fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
      fprintf(stderr,"      -U: Output sequence in upper case (default is lower case).\n");
      exit (1);
    }

  glen = strtod(argv[1],&eptr);
  if (*eptr != '\0')
    { fprintf(stderr,"%s: genome length is not a real number\n",Prog_Name);
      exit (1);
    }
  if (glen < 0.)
    { fprintf(stderr,"%s: Genome length must be positive (%g)\n",Prog_Name,glen);
      exit (1);
    }
  GENOME = (int) (glen*1000000.);

  //  Set up random number generator

  if (HASR)
    srand48(SEED);
  else
    srand48(getpid());

  //  Generate the sequence line at a time where all lines have width WDITH, save the last.


  fprintf(stdout,">random len=%d bias=%g\n",GENOME,BIAS);
  for (j = 0; j+WIDTH < GENOME; j += WIDTH)
    fprintf(stdout,"%s\n",random_genome(WIDTH));
  if (j < GENOME)
    fprintf(stdout,"%s\n",random_genome(GENOME-j));

  exit (0);
}
