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

#include "utilities.h"
#include "DB.h"

static char *Usage = "<path:db> ... ";

void HANDLER(char *path, char *name)
{ (void) name;
  unlink(path);
}

int main(int argc, char *argv[])
{ int   i;

  Prog_Name = Strdup("DBrm","");

  if (argc <= 1)
    fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);

  for (i = 1; i < argc; i++)
    if (List_DB_Files(argv[i],HANDLER))
      fprintf(stderr,"%s: Could not list database %s\n",Prog_Name,argv[i]);

  exit (0);
}
