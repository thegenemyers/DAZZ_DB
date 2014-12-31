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
 *  Interim code: upgrade previous db to have int's for pulse positions.
 *
 *  Author:  Gene Myers
 *  Date  :  September 2014
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

typedef struct
  { int     origin; //  Well #
    uint16  beg;    //  First pulse
    uint16  end;    //  Last pulse
    int64   boff;   //  Offset (in bytes) of compressed read in 'bases' file, or offset of
                    //    uncompressed bases in memory block
    int64   coff;   //  Offset (in bytes) of compressed quiva streams in 'quiva' file
    int     flags;  //  QV of read + flags above
  } HITS_OLD;

typedef struct
  { int     origin; //  Well #
    int     beg;    //  First pulse
    int     end;    //  Last pulse
    int64   boff;   //  Offset (in bytes) of compressed read in 'bases' file, or offset of
                    //    uncompressed bases in memory block
    int64   coff;   //  Offset (in bytes) of compressed quiva streams in 'quiva' file
    int     flags;  //  QV of read + flags above
  } HITS_NEW;

int main(int argc, char *argv[])
{ HITS_DB    db;
  FILE      *nxfile, *ixfile;
  char      *pwd, *root;
  int        i;

  if (argc != 2)
    { fprintf(stderr,"Usage: %s <path:db>\n",argv[0]);
      exit (1);
    }

  pwd    = PathTo(argv[1]);
  root   = Root(argv[1],".db");
  ixfile = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r");
  nxfile = Fopen(Catenate(pwd,PATHSEP,root,".ndx"),"w");
  if (ixfile == NULL || nxfile == NULL)
    exit (1);
  free(pwd);
  free(root);

  if (fread(&db,sizeof(HITS_DB),1,ixfile) != 1)
    SYSTEM_ERROR
  fwrite(&db,sizeof(HITS_DB),1,nxfile);

  for (i = 0; i < db.oreads; i++)
    { HITS_OLD  orec;
      HITS_NEW  nrec;

      if (fread(&orec,sizeof(HITS_OLD),1,ixfile) != 1)
        SYSTEM_ERROR

      nrec.origin = orec.origin;
      nrec.beg    = orec.beg;
      nrec.end    = orec.end;
      nrec.boff   = orec.boff;
      nrec.coff   = orec.coff;
      nrec.flags  = orec.flags;

      fwrite(&nrec,sizeof(HITS_NEW),1,nxfile);
    }

  fclose(ixfile);
  fclose(nxfile);

  exit (0);
}
