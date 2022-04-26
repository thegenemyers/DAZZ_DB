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
#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "DB.h"

static char *Usage = "[-vinf] <old:db|dam> <new:db|dam|dir>";

static int   VERBOSE;
static int   INQUIRE;
static int   NOXFER;
static int   FORCE;

static char *op;
static char  com[10];
static char *verb;

static char *nroot;
static char *npath;

  //  We assume this program uses so little code that memory allocation checks are unecessary

static void HANDLER(char *path, char *exten)
{ static char *cat = NULL;
  static int   max = -1;
  int   len;
  char *r;

  len = strlen(path) + strlen(npath) + strlen(nroot) + strlen(exten) + 50;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      cat = (char *) realloc(cat,max+1);
    }

  r = Root(path,"");
  if (*r == '.')
    { sprintf(cat,"%s %s %s/.%s.%s",com,path,npath,nroot,exten);
      if (VERBOSE)
        fprintf(stderr,"  %s %s to %s/.%s.%s\n",verb,path,npath,nroot,exten);
    }
  else
    { sprintf(cat,"%s %s %s/%s.%s",com,path,npath,nroot,exten);
      if (VERBOSE)
        fprintf(stderr,"  %s %s to %s/%s.%s\n",verb,path,npath,nroot,exten);
    }
  system(cat);
  free(r);
}

static int check(char *path, char *root)
{ int state;
  struct stat B;

  state = 0;
  if (stat(Catenate(path,"/",root,".db"),&B) == 0)
    state += 1;
  if (stat(Catenate(path,"/",root,".dam"),&B) == 0)
    state += 2;
  if (state == 3)
    { fprintf(stderr,"%s: %s refers to both a dam and a db ??\n",Prog_Name,root);
      exit (1);
    }
  return (state);
}

int main(int argc, char *argv[])
{ char *opath;
  char *oroot;
  int   otype,  ntype;
  struct stat B;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

#ifdef MOVE
    ARG_INIT("DBmv")
    op = "mv";
    verb = "Moving";
#else
    ARG_INIT("DBcp")
    op = "cp";
    verb = "Copying";
#endif

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vinf") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    INQUIRE = flags['i'];
    NOXFER  = flags['n'];
    FORCE   = flags['f'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  if (INQUIRE | NOXFER | FORCE)
    sprintf(com,"%s -%s%s%s",op,INQUIRE?"i":"",NOXFER?"n":"",FORCE?"f":"");
  else
    sprintf(com,"%s",op);

  opath = PathTo(argv[1]);
  if (strcmp(argv[1]+(strlen(argv[1])-3),".db") == 0)
    oroot = Root(argv[1],".db");
  else
    oroot = Root(argv[1],".dam");
  otype = check(opath,oroot);
  if (otype == 0)
    { fprintf(stderr,"%s: %s doesn't refer to a db or a dam\n",Prog_Name,oroot);
      exit (1);
    }

  if (stat(argv[2],&B) == 0 && (B.st_mode & S_IFMT) == S_IFDIR)
    { npath = argv[2];
      nroot = oroot;
      ntype = otype;
    }
  else
    { npath = PathTo(argv[2]);
      if (strcmp(argv[1]+(strlen(argv[1])-4),".db") == 0)
        nroot = Root(argv[2],".db");
      else
        nroot = Root(argv[2],".dam");
      ntype = check(npath,nroot);
      if (ntype != 0 && ntype != otype)
        { fprintf(stderr,"%s: one of %s and %s is a db, the other a dam ??\n",
                         Prog_Name,oroot,nroot);
          exit (1);
        }
    }

  if (List_DB_Files(argv[1],HANDLER) < 0)
    { fprintf(stderr,"%s: Could not find database %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  exit (0);
}
