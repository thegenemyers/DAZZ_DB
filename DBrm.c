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

static char *Usage = "[-vnf] <path:db|dam> ... ";

static int VERBOSE;
static int NODEL;
static int FORCE;

static char com[10];

static void HANDLER(char *path, char *exten)
{ static char *cat = NULL;
  static int   max = -1;
  int   len;

  (void) exten;

  len = strlen(path) + 50;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      cat = (char *) realloc(cat,max+1);
    }
  sprintf(cat,"%s %s",com,path);
  if (VERBOSE)
    fprintf(stderr,"  Deleting %s\n",path);
  system(cat);
}

int main(int argc, char *argv[])
{
  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("DBrm")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vnf") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    NODEL   = flags['n'];
    FORCE   = flags['f'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  if (NODEL | FORCE)
    sprintf(com,"rm -%s%s",NODEL?"n":"",FORCE?"f":"");
  else
    sprintf(com,"rm");

  { int i;

    for (i = 1; i < argc; i++)
      if (List_DB_Files(argv[i],HANDLER) < 0)
        fprintf(stderr,"%s: [WARNING] Could not find database %s\n",Prog_Name,argv[i]);
  }

  exit (0);
}
