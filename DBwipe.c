/*******************************************************************************************
 *
 *  Split a .db into a set of sub-database blocks for use by the Dazzler:
 *     Divide the database <path>.db conceptually into a series of blocks referable to on the
 *     command line as <path>.1.db, <path>.2.db, ...  If the -x option is set then all reads
 *     less than the given length are ignored, and if the -a option is not set then secondary
 *     reads from a given well are also ignored.  The remaining reads are split amongst the
 *     blocks so that each block is of size -s * 1Mbp except for the last which necessarily
 *     contains a smaller residual.  The default value for -s is 400Mbp because blocks of this
 *     size can be compared by our "overlapper" dalign in roughly 16Gb of memory.  The blocks
 *     are very space efficient in that their sub-index of the master .idx is computed on the
 *     fly when loaded, and the .bps file of base pairs is shared with the master DB.  Any
 *     tracks associated with the DB are also computed on the fly when loading a database block.
 *
 *  Author:  Gene Myers
 *  Date  :  September 2013
 *  Mod   :  New splitting definition to support incrementality, and new stub file format
 *  Date  :  April 2014
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

static char *Usage = "<path:db>";

int main(int argc, char *argv[])
{ HITS_DB    db;
  int        status;

  Prog_Name = Strdup("DBwipe","Allocating Program Name");

  if (argc != 2)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  //  Open db

  status = Open_DB(argv[1],&db);
  if (status < 0)
    exit (1);
  if (db.part > 0)
    { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  if (status)
    { fprintf(stderr,"%s: Cannot be called on a .dam: %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  { char    *pwd, *root;
    FILE    *index;
    int      i;

    pwd    = PathTo(argv[1]);
    root   = Root(argv[1],".db");
    unlink(Catenate(pwd,PATHSEP,root,".arw"));
    unlink(Catenate(pwd,PATHSEP,root,".qvs"));

    for (i = 0; i < db.nreads; i++)
      db.reads[i].coff = -1;
    db.allarr &= ~DB_ARROW;

    if ((index = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"w")) == NULL)
      { fprintf(stderr,"%s: Cannot open %s%s%s.idx\n",Prog_Name,pwd,PATHSEP,root);
        exit (1);
      }

    fwrite(&db,sizeof(HITS_DB),1,index);
    fwrite(db.reads,sizeof(HITS_READ),db.nreads,index);

    fclose(index);
  }

  Close_DB(&db);

  exit (0);
}
