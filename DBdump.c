/*******************************************************************************************
 *
 *  Display a portion of the data-base and selected information in 1-code format.
 *
 *  Author:  Gene Myers
 *  Date  :  November 2015
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage[] =
    { "[-rhsaiqf] [-uU] [-m<mask>]+",
      "           <path:db|dam> [ <reads:FILE> | <reads:range> ... ]"
    };

#define LAST_READ_SYMBOL   '$'
#define MAX_BUFFER       10001

typedef struct
  { FILE  *input;
    int    lineno;
    int    read;
    int    beg;
    int    end;
  } File_Iterator;

File_Iterator *init_file_iterator(FILE *input)
{ File_Iterator *it;

  it = Malloc(sizeof(File_Iterator),"Allocating file iterator");
  it->input = input;
  it->lineno = 1;
  rewind(input);
  return (it);
}

int next_read(File_Iterator *it)
{ static char nbuffer[MAX_BUFFER];

  char *eol;
  int   x;

  if (fgets(nbuffer,MAX_BUFFER,it->input) == NULL)
    { if (feof(it->input))
        return (1);
      SYSTEM_READ_ERROR;
    }
  if ((eol = index(nbuffer,'\n')) == NULL)
    { fprintf(stderr,"%s: Line %d in read list is longer than %d chars!\n",
                     Prog_Name,it->lineno,MAX_BUFFER-1);
      return (1);
    }
  *eol = '\0';
  x = sscanf(nbuffer," %d %d %d",&(it->read),&(it->beg),&(it->end));
  if (x == 1)
    it->beg = -1;
  else if (x != 3)
    { fprintf(stderr,"%s: Line %d of read list is improperly formatted\n",Prog_Name,it->lineno);
      return (1);
    }
  it->lineno += 1;
  return (0);
}

static int qv_map[51] =
  { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
    'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D',
    'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
    'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
    'Y'
  };

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  int         Quiva_DB, Arrow_DB;
  int         FirstRead;
  FILE       *hdrs      = NULL;
  char       *hdrs_name = NULL;

  int         nfiles;
  char      **fhead = NULL;
  char      **ffile = NULL;
  int        *findx = NULL;
  char       *empty = "";

  int            input_pts;
  int            reps = 0;
  int           *pts  = NULL;
  File_Iterator *iter = NULL;
  FILE          *input = NULL;

  int         TRIM, UPPER;
  int         DORED, DOSEQ, DOARW, DOQVS, DOHDR, DOIQV, DOFLN, DAM;

  int          MMAX, MTOP;
  char       **MASK;
  DAZZ_TRACK **MTRACK;
  int       **MDATA;

  DAZZ_TRACK *qv_track;
  uint8      *qv_data = NULL;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("DBdump")

    MTOP  = 0;
    MMAX  = 10;
    MASK  = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    if (MASK == NULL)
      exit (1);

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("hfqrsaiuU")
            break;
          case 'm':
            if (MTOP >= MMAX)
              { MMAX = 1.2*MTOP + 10;
                MASK = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
                if (MASK == NULL)
                  exit (1);
              }
            MASK[MTOP++] = argv[i]+2;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    DAM   = 0;
    TRIM  = 1-flags['u'];
    UPPER = 1+flags['U'];
    DOQVS = flags['q'];
    DORED = flags['r'];
    DOSEQ = flags['s'];
    DOARW = flags['a'];
    DOHDR = flags['h'];
    DOFLN = flags['f'];
    DOIQV = flags['i'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -r: R #              - read number\n");
        fprintf(stderr,"      -h: H # string       - fasta header prefix\n");
        fprintf(stderr,"          L # # #          - location: well, pulse start, pulse end\n");
        fprintf(stderr,"          Q #              - quality of read (#/1000)\n");
        fprintf(stderr,"      -s: S # string       - sequence string\n");
        fprintf(stderr,"      -a: N # # # #        - SNR of ACGT channels (#/100)\n");
        fprintf(stderr,"          A # string       - arrow pulse-width string\n");
        fprintf(stderr,"      -i: I # string       ");
        fprintf(stderr,"- intrinsic quality vector (as an ASCII string)\n");
        fprintf(stderr,"      -q: d # string       - Quiva deletion values (as an ASCII string)\n");
        fprintf(stderr,"          c # string       - Quiva deletion character string\n");
        fprintf(stderr,"          i # string       - Quiva insertion value string\n");
        fprintf(stderr,"          m # string       - Quiva merge value string\n");
        fprintf(stderr,"          s # string       - Quiva substitution value string\n");
        fprintf(stderr,"      -f: F # string       - file name for all until next F\n");
        fprintf(stderr,"      -m: Tx #n (#b #e)^#n ");
        fprintf(stderr,"- x'th track on command line, #n intervals all on same line\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -u: Dump entire untrimmed database.\n");
        fprintf(stderr,"      -U: Output base pairs in upper case letters\n");
        exit (1);
      }

    if ( ! TRIM && DOIQV)
      { fprintf(stderr,"%s: -i and -u are incompatible\n",Prog_Name);
        exit (1);
      }
  }

  //  Open DB or DAM, and if a DAM open also .hdr file

  { char *pwd, *root;
    int   status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 1)
      { root   = Root(argv[1],".dam");
        pwd    = PathTo(argv[1]);

        if (db->part > 0)
          *rindex(root,'.') = '\0';
        hdrs_name = Strdup(Catenate(pwd,PATHSEP,root,".hdr"),"Allocating header file name");
        hdrs      = Fopen(hdrs_name,"r");
        if (hdrs == NULL)
          exit (1);
        DAM = 1;
        if (DOQVS)
          { fprintf(stderr,"%s: -q option is not compatible with a .dam DB\n",Prog_Name);
            exit (1);
          }
        if (DOARW)
          { fprintf(stderr,"%s: -a option is not compatible with a .dam DB\n",Prog_Name);
            exit (1);
          }

        free(root);
        free(pwd);
      }

    Arrow_DB = ((db->allarr & DB_ARROW) != 0);
    Quiva_DB = (db->reads[0].coff >= 0 && (db->allarr & DB_ARROW) == 0);
    if (DOARW)
      { if (!Arrow_DB)
          { fprintf(stderr,"%s: -a option set but no Arrow data in DB\n",Prog_Name);
            exit (1);
          }
      }
    if (DOQVS)
      { if (!Quiva_DB)
          { fprintf(stderr,"%s: -q option set but no Quiver data in DB\n",Prog_Name);
            exit (1);
          }
      }
  }

  //  Load QVs if requested

  if (DOQVS)
    { if (Open_QVs(db) < 0)
        { fprintf(stderr,"%s: QVs requested, but no .qvs for data base\n",Prog_Name);
          exit (1);
	}
    }
  if (DOARW)
    { if (Open_Arrow(db) < 0)
        { fprintf(stderr,"%s: Arrow requested, but no .arr for data base\n",Prog_Name);
          exit (1);
	}
    }

  //  Check tracks and load tracks for untrimmed DB

  { int i, status, kind;

    MTRACK = Malloc(sizeof(DAZZ_TRACK *)*MTOP,"Allocation of track pointer vector");
    MDATA  = Malloc(sizeof(int *)*MTOP,"Allocation of track data buffer pointers");
    if (MTRACK == NULL || MDATA == NULL)
      exit (1);

    for (i = 0; i < MTOP; i++)
      { status = Check_Track(db,MASK[i],&kind);
        if (status == -2)
          { fprintf(stderr,"%s: Warning: -m%s option given but no track found.\n",
                           Prog_Name,MASK[i]);
            exit (1);
          }
        else if (status == -1)
          { fprintf(stderr,"%s: Warning: %s track not sync'd with db.\n",Prog_Name,MASK[i]);
            exit (1);
          }
        else if (kind != MASK_TRACK)
          { fprintf(stderr,"%s: Warning: %s track is not a mask track.\n",Prog_Name,MASK[i]);
            exit (1);
          }
        else if (status == 0)
          MTRACK[i] = Open_Track(db,MASK[i]);
        else if (status == 1 && !TRIM)
          { fprintf(stderr,"%s: Warning: %s track is for a trimmed db but -u is set.\n",
                           Prog_Name,MASK[i]);
            exit (1);
          }
      }
  }

  //  If get prolog and file names and index ranges from the .db or .dam file 

  { char *pwd, *root;
    FILE *dstub;
    char *dstub_name;
    int   i;

    if (DAM)
      root = Root(argv[1],".dam");
    else
      root = Root(argv[1],".db");
    pwd = PathTo(argv[1]);
    if (db->part > 0)
      *rindex(root,'.') = '\0';
    if (DAM)
      dstub_name = Strdup(Catenate(pwd,"/",root,".dam"),"Allocating dam file name");
    else
      dstub_name = Strdup(Catenate(pwd,"/",root,".db"),"Allocating db file name");
    dstub = Fopen(dstub_name,"r");
    if (dstub_name == NULL || dstub == NULL)
      exit (1);
    free(pwd);
    free(root);

    FSCANF(dstub,DB_NFILE,&nfiles)

    fhead = (char **) Malloc(sizeof(char *)*nfiles,"Allocating file list");
    ffile = (char **) Malloc(sizeof(char *)*(nfiles+1),"Allocating file list");
    findx = (int *) Malloc(sizeof(int *)*(nfiles+1),"Allocating file index");
    if (fhead == NULL || ffile == NULL || findx == NULL)
      exit (1);

    findx += 1;
    findx[-1] = 0;
    ffile += 1;
    ffile[-1] = empty;

    for (i = 0; i < nfiles; i++)
      { char prolog[MAX_NAME], fname[MAX_NAME];

        FSCANF(dstub,DB_FDATA,findx+i,fname,prolog)
        if ((fhead[i] = Strdup(prolog,"Adding to file list")) == NULL)
          exit (1);
        if ((ffile[i] = Strdup(fname,"Adding to file list")) == NULL)
          exit (1);
      }

    free(dstub_name);
    fclose(dstub);

    //  If TRIM (the default) then "trim" prolog ranges and the DB

    if (TRIM)
      { int        nid, oid, lid;
        int        cutoff, allflag;
        DAZZ_READ *reads;

        reads  = db->reads - db->ufirst;
        cutoff = db->cutoff;
        if ((db->allarr & DB_ALL) != 0)
          allflag = 0;
        else
          allflag = DB_BEST;
        
        nid = 0;
        oid = db->ufirst;
        lid = oid + db->nreads;
        for (i = 0; i < nfiles; i++)
          { while (oid < findx[i] && oid < lid)
              { if ((reads[oid].flags & DB_BEST) >= allflag && reads[oid].rlen >= cutoff)
                  nid++;
                oid += 1;
              }
            findx[i] = nid;
          }
      }

    else if (db->part > 0)
      { for (i = 0; i < nfiles; i++)
          findx[i] -= db->ufirst;
      }
  }

  if (TRIM)
    { int i, status, kind;

      Trim_DB(db);

      //  Load tracks for trimmed DB

      for (i = 0; i < MTOP; i++)
        { status = Check_Track(db,MASK[i],&kind);
          if (status > 0)
            MTRACK[i] = Open_Track(db,MASK[i]);
        }
      FirstRead = db->tfirst;
    }
  else
    FirstRead = db->ufirst;

  { int c;

    for (c = 0; c < MTOP; c++)
      MDATA[c] = (int *) New_Track_Buffer(MTRACK[c]);
  }

  if (DOIQV)
    { int status, kind;

      status = Check_Track(db,"qual",&kind);
      if (status == -2)
        { fprintf(stderr,"%s: .qual-track does not exist for this db.\n",Prog_Name);
          exit (1);
        }
      if (status == -1)
        { fprintf(stderr,"%s: .qual-track not sync'd with db.\n",Prog_Name);
          exit (1);
        }
      qv_track = Open_Track(db,"qual");
      qv_data  = (uint8 *) New_Track_Buffer(qv_track);
    }

  //  Process read index arguments into a list of read ranges

  input_pts = 0;
  if (argc == 3)
    { if (argv[2][0] != LAST_READ_SYMBOL || argv[2][1] != '\0')
        { char *eptr, *fptr;
          int   b, e;

          b = strtol(argv[2],&eptr,10);
          if (eptr > argv[2] && b > 0)
            { if (*eptr == '-')
                { if (eptr[1] != LAST_READ_SYMBOL || eptr[2] != '\0')
                    { e = strtol(eptr+1,&fptr,10);
                      input_pts = (fptr <= eptr+1 || *fptr != '\0' || e <= 0);
                    }
                }
              else
                input_pts = (*eptr != '\0');
            }
          else
            input_pts = 1;
        }
    }

  if (input_pts)
    { input = Fopen(argv[2],"r");
      if (input == NULL)
        exit (1);

      iter = init_file_iterator(input);
    }
  else
    { pts = (int *) Malloc(sizeof(int)*2*(argc-1),"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      reps = 0;
      if (argc > 2)
        { int   c, b, e;
          char *eptr, *fptr;

          for (c = 2; c < argc; c++)
            { if (argv[c][0] == LAST_READ_SYMBOL)
                { b = db->nreads;
                  eptr = argv[c]+1;
                }
              else
                b = strtol(argv[c],&eptr,10);
              if (eptr > argv[c])
                { if (b <= 0)
                    { fprintf(stderr,"%s: %d is not a valid index\n",Prog_Name,b);
                      exit (1);
                    }
                  if (*eptr == 0)
                    { pts[reps++] = b;
                      pts[reps++] = b;
                      continue;
                    }
                  else if (*eptr == '-')
                    { if (eptr[1] == LAST_READ_SYMBOL)
                        { e = db->nreads;
                          fptr = eptr+2;
                        }
                      else
                        e = strtol(eptr+1,&fptr,10);
                      if (fptr > eptr+1 && *fptr == 0 && e > 0)
                        { pts[reps++] = b;
                          pts[reps++] = e;
                          if (b > e)
                            { fprintf(stderr,"%s: Empty range '%s'\n",Prog_Name,argv[c]);
                              exit (1);
                            }
                          continue;
                        }
                    }
                }
              fprintf(stderr,"%s: argument '%s' is not an integer range\n",Prog_Name,argv[c]);
              exit (1);
            }
        }
      else
        { pts[reps++] = 1;
          pts[reps++] = db->nreads;
        }
    }

  //  Scan to count the size of things

  { DAZZ_READ  *reads;
    int         c, b, e, i, m;
    int         map, substr, last;
    int64       noreads;
    int64       seqmax, seqtot;
    int64       iqvmax, iqvtot;
    int64       flnmax, flntot;
    int64       hdrmax, hdrtot;
    int64       trkmax[MTOP], trktot[MTOP];

    map    = 0;
    reads  = db->reads;
    substr = 0;
    last   = -1;

    noreads = 0;
    seqmax = 0;
    seqtot = 0;
    iqvmax = 0;
    iqvtot = 0;
    flnmax = 0;
    flntot = 0;
    hdrmax = 0;
    hdrmax = 0;
    hdrtot = 0;
    for (m = 0; m < MTOP; m++)
      { trkmax[m] = 0;
        trktot[m] = 0;
      }

    c = 0;
    while (1)
      { if (input_pts)
          { if (next_read(iter))
              break;
            e = iter->read;
            b = e-1;
            substr = (iter->beg >= 0);
          }
        else
          { if (c >= reps)
              break;
            b = pts[c]-1;
            e = pts[c+1];
            if (e > db->nreads)
              e = db->nreads;
            c += 2;
          }
 
        if (map > 0 && findx[map-1] <= b+FirstRead && b+FirstRead < findx[map])
          ;
        else
          { map = 0;
            while (b + FirstRead >= findx[map])
              map += 1;
            map -= 1;
          }

        for (i = b; i < e; i++)
          { int         len, ten;
            int         fst, lst;
            DAZZ_READ  *r;

            r   = reads + i;
            len = r->rlen;

            noreads += 1;

            if (DOFLN && i+FirstRead >= findx[map])
              { int ten;

                if (strcmp(ffile[map+1],ffile[last]) != 0)
                  { ten = strlen(ffile[map+1]);
                    if (flnmax < ten)
                      flnmax = ten;
                    flntot += ten;
                    last = map+1;
                  }
                if (!DOHDR || DAM)
                  map += 1;
              }

            if (DOHDR)
              { int ten;

                if (DAM)
                  { char header[MAX_NAME];

                    FSEEKO(hdrs,r->coff,SEEK_SET)
                    FGETS(header,MAX_NAME,hdrs)
                    header[strlen(header)-1] = '\0';
                    ten = strlen(header);

                    if (hdrmax < ten)
                      hdrmax = ten;
                    hdrtot += ten;
                  }
                else if (i+FirstRead >= findx[map])
                  { map += 1;
                    ten = strlen(fhead[map]);

                    if (hdrmax < ten)
                      hdrmax = ten;
                    hdrtot += ten;
                  }
              }

            for (m = 0; m < MTOP; m++)
              { ten = MTRACK[m]->alen[i];
                if (ten > trkmax[m])
                  trkmax[m] = ten;
                trktot[m] += ten;
              }

            if (substr)
              { fst = iter->beg;
                lst = iter->end;
                if (DOIQV)
                  { fprintf(stderr,"%s: Cannot select subreads when -i is requested\n",Prog_Name);
                    exit (1);
                  }
              }
            else
              { fst = 0;
                lst = len;
              }

            if (DOSEQ | DOQVS | DOARW)
              { int ten = lst-fst;
                if (ten > seqmax)
                  seqmax = ten;
                seqtot += ten;
              }
            if (DOIQV)
              { int ten = qv_track->alen[i];
                if (ten > iqvmax)
                  iqvmax = ten;
                iqvtot += ten;
              }
          }
      }

    PRINTF("+ R %lld\n",noreads)
    PRINTF("+ M %d\n",MTOP)
    if (DOHDR)
      { PRINTF("+ H %lld\n",hdrtot)
        PRINTF("@ H %lld\n",hdrmax)
      }
    for (m = 0; m < MTOP; m++)
      { PRINTF("+ T%d %lld\n",m,trktot[m])
        PRINTF("@ T%d %lld %ld %s\n",m,trkmax[m],strlen(MASK[m]),MASK[m])
      }
    if (DOSEQ | DOQVS | DOARW)
      { PRINTF("+ S %lld\n",seqtot)
        PRINTF("@ S %lld\n",seqmax)
      }
    if (DOIQV)
      { PRINTF("+ I %lld\n",iqvtot)
        PRINTF("@ I %lld\n",iqvmax)
      }
    if (DOFLN)
      { PRINTF("+ F %lld\n",flntot)
        PRINTF("@ F %lld\n",flnmax)
      }
  }

  //  Display each read (and/or QV streams) in the active DB according to the
  //    range pairs in pts[0..reps) and according to the display options.

  { DAZZ_READ  *reads;
    char       *read, *arrow, **entry;
    int         c, b, e, i, m;
    int         substr;
    int         map, last;
    char        qvname[5] = { 'd', 'c', 'i', 'm', 's' };

    read  = New_Read_Buffer(db);
    if (DOQVS)
      entry = New_QV_Buffer(db);
    else
      entry = NULL;
    if (DOARW)
      arrow = New_Read_Buffer(db);
    else
      arrow = NULL;

    map    = 0;
    reads  = db->reads;
    substr = 0;
    last   = -1;

    if (input_pts)
      iter = init_file_iterator(input);
    else
      iter = NULL;

    c = 0;
    while (1)
      { if (input_pts)
          { if (next_read(iter))
              break;
            e = iter->read;
            b = e-1;
            substr = (iter->beg >= 0);
          }
        else
          { if (c >= reps)
              break;
            b = pts[c]-1;
            e = pts[c+1];
            if (e > db->nreads)
              e = db->nreads;
            c += 2;
          }

        if (map > 0 && findx[map-1] <= b+FirstRead && b+FirstRead < findx[map])
          ;
        else
          { map = 0;
            while (b + FirstRead >= findx[map])
              map += 1;
            map -= 1;
          }

        for (i = b; i < e; i++)
          { int         len;
            int         fst, lst;
            int         flags, qv;
            DAZZ_READ  *r;

            r   = reads + i;
            len = r->rlen;
            if (DORED)
              printf("R %d\n",FirstRead + (i+1));

            flags = r->flags;
            qv    = (flags & DB_QV);

            if (DOFLN && i+FirstRead >= findx[map])
              { if (strcmp(ffile[map+1],ffile[last]) != 0)
                  { PRINTF("F %ld %s\n",strlen(ffile[map+1]),ffile[map+1])
                    last = map+1;
                  }
                if (!DOHDR || DAM)
                  map += 1;
              }

            if (DOHDR)
              { if (DAM)
                  { char header[MAX_NAME];

                    FSEEKO(hdrs,r->coff,SEEK_SET)
                    FGETS(header,MAX_NAME,hdrs)
                    header[strlen(header)-1] = '\0';
                    PRINTF("H %ld %s\n",strlen(header),header)
                    PRINTF("L %d %d %d\n",r->origin,r->fpulse,r->fpulse+len)
                  }
                else
                  { if (i+FirstRead >= findx[map])
                      { map += 1;
                        PRINTF("H %ld %s\n",strlen(fhead[map]),fhead[map])
                      }
                    PRINTF("L %d %d %d\n",r->origin,r->fpulse,r->fpulse+len)
                    if (Quiva_DB && qv > 0)
                      PRINTF("Q %d\n",qv)
                    else if (Arrow_DB)
                      { int   j, snr[4];
                        int64 big;

                        big = *((uint64 *) &(r->coff));
                        for (j = 0; j < 4; j++)
                          { snr[3-j] = (big & 0xffff);
                            big    >>= 16;
                          }
                        PRINTF("N %d %d %d %d\n",snr[0],snr[1],snr[2],snr[3])
                      }
                  }
              }

            if (DOQVS)
              Load_QVentry(db,i,entry,UPPER);
            if (DOSEQ)
              Load_Read(db,i,read,UPPER);
            if (DOARW)
              Load_Arrow(db,i,arrow,1);

            for (m = 0; m < MTOP; m++)
              { int   *d;
                int64  f, j;

                d = MDATA[m];
                f = (Load_Track_Data(MTRACK[m],i,d) >> 2);

                PRINTF("T%d %lld",m,f/2)
                if (f > 0)
                  { for (j = 0; j < f; j += 2)
                      PRINTF(" %d %d",d[j],d[j+1])
                  }
                PRINTF("\n")
              }

            if (substr)
              { fst = iter->beg;
                lst = iter->end;
              }
            else
              { fst = 0;
                lst = len;
              }

            if (DOSEQ)
              { PRINTF("S %d ",lst-fst)
                PRINTF("%.*s\n",lst-fst,read+fst)
              }

            if (DOARW)
              { PRINTF("A %d ",lst-fst)
                PRINTF("%.*s\n",lst-fst,arrow+fst)
              }

            if (DOIQV)
              { int j, f;

                f = Load_Track_Data(qv_track,i,qv_data);
                PRINTF("I %d ",f)
                for (j = 0; j < f; j++)
                  { if (putchar(qv_map[qv_data[j]]) == EOF)
                      SYSTEM_WRITE_ERROR
                  }
                PRINTF("\n")
              }

            if (DOQVS)
              { int k;

                for (k = 0; k < 5; k++)
                  { PRINTF("%c %d ",qvname[k],lst-fst)
                    PRINTF("%.*s\n",lst-fst,entry[k]+fst)
                  }
              }
          }
      }
  }

  FCLOSE(stdout)

  if (input_pts)
    { fclose(input);
      free(iter);
    }
  else
    free(pts);

  if (DAM)
    fclose(hdrs);
  else
    { int i;

      for (i = 0; i < nfiles; i++)
        { free(fhead[i]);
          free(ffile[i]);
        }
      free(fhead);
      free(ffile-1);
      free(findx-1);
    }
  Close_DB(db);

  exit (0);
}
