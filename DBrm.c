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

static char *Usage = "[-v] <path:db|dam> ... ";

static int VERBOSE;

static void HANDLER(char *path, char *exten)
{ (void) exten;
  if (unlink(path) != 0)
    fprintf(stderr,"%s: [WARNING] Couldn't delete file %s\n",Prog_Name,path);
  else if (VERBOSE)
    fprintf(stderr,"  Deleting %s\n",path);
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
        { ARG_FLAGS("v") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { int i;

    for (i = 1; i < argc; i++)
      if (List_DB_Files(argv[i],HANDLER) < 0)
        fprintf(stderr,"%s: [WARNING] Could not find database %s\n",Prog_Name,argv[i]);
  }

  exit (0);
}
