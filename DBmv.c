/********************************************************************************************
 *
 *  Remove a list of .db databases
 *     Delete all the files for the given data bases <name>.db ... (there are a couple
 *     of hidden . files for each DB, and these are removed too.)  Do not use "rm" to
 *     remove a database.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "DB.h"

static char *Usage = "[-v] <old:db|dam> <new:db|dam>";

static int   VERBOSE;
static char *nroot;
static char *npath;

  //  We assume this program uses so little code that memory allocation checks are unecessary

static char *Catenate5(char *path, char *sep1, char *root, char *sep2, char *suffix)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  len = strlen(path) + strlen(sep1) + strlen(root) + strlen(sep2) + strlen(suffix);
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      cat = (char *) realloc(cat,max+1);
    }
  sprintf(cat,"%s%s%s%s%s",path,sep1,root,sep2,suffix);
  return (cat);
}

static void HANDLER(char *path, char *exten)
{ char *r, *n;
  
  r = Root(path,"");
  if (*r == '.')
    n = Catenate5(npath,"/.",nroot,".",exten);
  else
    n = Catenate5(npath,"/",nroot,".",exten);
  if (rename(path,n) != 0)
    fprintf(stderr,"%s: [WARNING] Couldn't rename file %s\n",Prog_Name,r);
  else if (VERBOSE)
    fprintf(stderr,"  Moving %s to %s\n",path,n);
  free(r);
}

int main(int argc, char *argv[])
{
  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("DBmv")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("v") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  if (strcmp(argv[1]+(strlen(argv[1])-4),".dam") == 0)
    nroot = Root(argv[2],".dam");
  else
    nroot = Root(argv[2],".db");
  npath = PathTo(argv[2]);

printf(" From = '%s'\n",argv[1]);
  if (List_DB_Files(argv[1],HANDLER) < 0)
    { fprintf(stderr,"%s: Could not find database %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  exit (0);
}
