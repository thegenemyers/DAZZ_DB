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

typedef struct
  { int     origin; //  Well #
    int     beg;    //  First pulse
    int     end;    //  Last pulse
    int64   boff;   //  Offset (in bytes) of compressed read in 'bases' file, or offset of
                    //    uncompressed bases in memory block
    int64   coff;   //  Offset (in bytes) of compressed quiva streams in 'quiva' file
    int     flags;  //  QV of read + flags above
  } HITS_OLD;

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

  fprintf(stderr,"Converting %d reads\n",db.ureads);
  fflush(stderr);

  for (i = 0; i < db.ureads; i++)
    { HITS_OLD  orec;
      HITS_READ nrec;

      if (i%10000 == 0)
        { fprintf(stderr,"  Processing %d\n",i);
          fflush(stderr);
        }

      if (fread(&orec,sizeof(HITS_OLD),1,ixfile) != 1)
        SYSTEM_ERROR

      nrec.origin = orec.origin;
      nrec.fpulse = orec.beg;
      nrec.rlen   = orec.end-orec.beg;
      nrec.boff   = orec.boff;
      nrec.coff   = orec.coff;
      nrec.flags  = orec.flags;

      fwrite(&nrec,sizeof(HITS_READ),1,nxfile);
    }

  fclose(ixfile);
  fclose(nxfile);

  exit (0);
}
