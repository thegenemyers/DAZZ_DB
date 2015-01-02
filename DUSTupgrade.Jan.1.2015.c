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
 *  Interim code: upgrade previous db to have fpulse,rlen fields
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
