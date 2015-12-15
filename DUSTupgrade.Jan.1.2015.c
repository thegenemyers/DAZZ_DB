/*******************************************************************************************
 *
 *  Interim code: upgrade dust track indices from int's to int64's
 *
 *  Author:  Gene Myers
 *  Date  :  December 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

int main(int argc, char *argv[])
{ FILE      *afile, *dfile;
  FILE      *nafile, *ndfile;
  char      *pwd, *root;
  int        size, tracklen;
  int        i, vint, dint;
  int64      vlong;

  if (argc != 2)
    { fprintf(stderr,"Usage: %s <path:db>\n",argv[0]);
      exit (1);
    }

  pwd    = PathTo(argv[1]);
  root   = Root(argv[1],".db");
  afile  = Fopen(Catenate(pwd,PATHSEP,root,".dust.anno"),"r");
  dfile  = Fopen(Catenate(pwd,PATHSEP,root,".dust.data"),"r");
  nafile = Fopen(Catenate(pwd,PATHSEP,root,".next.anno"),"w");
  ndfile = Fopen(Catenate(pwd,PATHSEP,root,".next.data"),"w");
  if (afile == NULL || dfile == NULL || nafile == NULL || ndfile == NULL)
    exit (1);
  free(pwd);
  free(root);

  if (fread(&tracklen,sizeof(int),1,afile) != 1)
    SYSTEM_ERROR
  fwrite(&tracklen,sizeof(int),1,nafile);

  if (fread(&size,sizeof(int),1,afile) != 1)
    SYSTEM_ERROR
  size = 8;
  fwrite(&size,sizeof(int),1,nafile);

  for (i = 0; i <= tracklen; i++)
    { if (fread(&vint,sizeof(int),1,afile) != 1)
        SYSTEM_ERROR
      vlong = vint;
      fwrite(&vlong,sizeof(int64),1,nafile);
    }

  vint >>= 2;
  for (i = 0; i < vint; i += 2)
    { if (fread(&dint,sizeof(int),1,dfile) != 1)
        SYSTEM_ERROR
      fwrite(&dint,sizeof(int),1,ndfile);
      if (fread(&dint,sizeof(int),1,dfile) != 1)
        SYSTEM_ERROR
      dint += 1;
      fwrite(&dint,sizeof(int),1,ndfile);
    }

  fclose(nafile);
  fclose(ndfile);
  fclose(afile);
  fclose(dfile);

  exit (0);
}
