#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "DB.h"

static char *Usage = "-v [-s<int(400)>] [-c<int(6000)>] <path:db>";

int main(int argc, char *argv[]) { 
    HITS_DB    _db, *db = &_db;
    int         VERBOSE, SIZE, CUTOFF;
    FILE       *dbfile;
    int         nfiles;

    //  Process arguments

    { 
        int   i, j;
        int   flags[128];
        char *eptr;

        ARG_INIT("DB2idmap")

        SIZE   = 400;
        CUTOFF = 6000;
        j = 1;
        for (i = 1; i < argc; i++) {
            if (argv[i][0] == '-')
                switch (argv[i][1])
                { default:
                    break;
                    case 's':
                        ARG_POSITIVE(SIZE,"Block size")
                        break;
                    case 'c':
                        ARG_POSITIVE(CUTOFF,"Seed length cutoff")
                        break;
                }
            else
                argv[j++] = argv[i];
        }
        argc = j;

        VERBOSE = flags['v'];

        if (argc != 2) { 
            fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
            exit (1);
        }
    }

    //  Open db and also db image file (dbfile)

    if (Open_DB(argv[1],db)) { 
        fprintf(stderr,"%s: Database %s.db could not be opened\n",Prog_Name,argv[1]);
        exit (1);
    }

    if (VERBOSE) {
        fprintf(stderr,"Obtaining read id mapping information ...\n");
    }

    char  *basename;
    { char *pwd;

        pwd    = PathTo(argv[1]);
        basename = Root(argv[1],".db");
        dbfile = Fopen(Catenate(pwd,"/",basename,".db"),"r");
        free(pwd);
        if (dbfile == NULL)
            exit (1);

    }

    if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
        SYSTEM_ERROR

    { 
        HITS_READ  *reads;
        int         i, fcount, last, nblock, ireads, breads;
        int64       size, totlen;
        char  prolog[MAX_NAME], fname[MAX_NAME];

        size = SIZE*1000000ll;

        reads = db->reads;
        
        ireads = 0;
        breads = 0;
        totlen = 0;
        nblock = 1;
        fcount = 0;

        if (fscanf(dbfile,DB_FDATA,&last,fname,prolog) != 3)
          SYSTEM_ERROR

        FILE *ofile;
        if ((ofile = Fopen(Catenate(".","/",basename,".idmap"),"w")) == NULL)
            exit (1);

        for (i = 0; i < db->nreads; i++)
        {   
            int        len, flags;
            HITS_READ *r;

            r     = reads + i;
            len   = r->end - r->beg;
            flags = r->flags;

            // same logic from DBsplit.c
            if (len >= db->cutoff && (flags & DB_BEST) != 0) {
                ireads += 1;
                breads += 1;
                totlen += len;
    
                if (totlen >= size || ireads >= READMAX) {
                    ireads = 0;
                    totlen = 0;
                    nblock += 1;
                }
                unsigned short seed = len>=CUTOFF?1:0;
                fprintf(ofile,"%d %d %s/%d/%d_%d %d\n",
                        nblock,breads,prolog,r->origin,r->beg,r->end,seed);
            }

            if (i+1 >= last && ++fcount < nfiles) {
                if (fscanf(dbfile,DB_FDATA,&last,fname,prolog) != 3)
                    SYSTEM_ERROR
            }

        }
        fclose(ofile);
    }
    fclose(dbfile);
    Close_DB(db);

    exit (0);
}
