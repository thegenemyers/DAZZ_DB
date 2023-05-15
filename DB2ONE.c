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
#include "ONElib.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage[] =
    { "[-u] [-aqhwf] [-m<mask>]+",
      "           <path:db|dam> [ <reads:FILE> | <reads:range> ... ]"
    };

static char *One_Schema =
  "P 3 daz                   This is a Dazzler DB 1-code file\n"

  "G f 2 3 INT 6 STRING      Source file reads came from\n"

  "O R 2 3 INT 6 STRING      Read idx and string for read\n"

  "D A 1 6 STRING                 Arrow pulse widths\n"

  "D D 1 6 STRING                 Quiva del vals\n"
  "D C 1 6 STRING                 Quiva del char\n"
  "D I 1 6 STRING                 Quiva ins vals\n"
  "D M 1 6 STRING                 Quiva mrg vals\n"
  "D S 1 6 STRING                 Quiva sub vals\n"

  "D H 1 6 STRING                 Original fasta/q header\n"

  "D W 3 3 INT 3 INT 3 INT           well, pulse start, pulse end (for db's)\n"
  "D N 4 3 INT 3 INT 3 INT 3 INT        SNR of ACGT channels (if Arrow-DB)\n"
  "D Q 1 3 INT                          read quality (if Quiva-DB)\n"

  "D G 3 3 INT 3 INT 3 INT           contig, firstbp, lastbp (for dam's)\n"

  "D X 2 3 INT 6 STRING            Prolog: name of track idx\n"
  "D T 2 3 INT 8 INT_LIST          Track idx, interval pairs list\n";


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

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  DAZZ_STUB  *stub;
  int         Quiva_DB, Arrow_DB;
  int         FirstRead;
  FILE       *hdrs      = NULL;
  char       *hdrs_name = NULL;

  OneSchema *schema;
  OneFile   *file1;
  char      *command;

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

  int         TRIM, DAM;
  int         DOARW, DOQVS, DOHDR, DOFLN, DOWEL;

  int          MMAX, MTOP;
  char       **MASK;
  DAZZ_TRACK **MTRACK;
  int       **MDATA;
  int64      *MBUFFER;

  //  Process arguments and capture command line for provenance

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("DB2ONE")

    { int   n, t;
      char *c;

      n = 0;
      for (t = 1; t < argc; t++)
        n += strlen(argv[t])+1;

      command = Malloc(n+1,"Allocating command string");
      if (command == NULL)
        exit (1);

      c = command;
      if (argc >= 1)
        { c += sprintf(c,"%s",argv[1]);
          for (t = 2; t < argc; t++)
            c += sprintf(c," %s",argv[t]);
        }
      *c = '\0';
    }

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
            ARG_FLAGS("uaqhwf")
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
    DOARW = flags['a'];
    DOQVS = flags['q'];
    DOHDR = flags['h'];
    DOWEL = flags['w'];
    DOFLN = flags['f'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      Output read idx and sequence by default (R line)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -a: Output truncated arrow pulse-width string (A line)\n");
        fprintf(stderr,"      -q: Quiver edit vectors (D, C, I, M, and S lines)\n");
        fprintf(stderr,"      -h: Output fasta header prefix (H line)\n");
        fprintf(stderr,"      -w: Output well, pulse start and end (if .db, W line)\n");
        fprintf(stderr,"            + SNR of ACGT channels (if Arrow DB, N line)\n");
        fprintf(stderr,"            + quality value of read (if Quiver DB, Q line)\n");
        fprintf(stderr,"      -w: Contig, firstbp, and lastbp (if .dam, G line)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -f: group by origin file (f line)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -m: output *mask* track (T line)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -u: Dump entire untrimmed database.\n");
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
    int   i;

    if (DAM)
      root = Root(argv[1],".dam");
    else
      root = Root(argv[1],".db");
    pwd = PathTo(argv[1]);
    if (db->part > 0)
      *rindex(root,'.') = '\0';
    if (DAM)
      stub = Read_DB_Stub(Catenate(pwd,"/",root,".dam"),
                           DB_STUB_NREADS|DB_STUB_FILES|DB_STUB_PROLOGS);
    else
      stub = Read_DB_Stub(Catenate(pwd,"/",root,".db"),
                           DB_STUB_NREADS|DB_STUB_FILES|DB_STUB_PROLOGS);
    free(pwd);
    free(root);

    fhead  = stub->prolog;
    ffile  = stub->fname;
    findx  = stub->nreads;
    nfiles = stub->nfiles;

    findx[-1] = 0;
    ffile[-1] = empty;

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

  { int c, n;

    n = 0;
    for (c = 0; c < MTOP; c++)
      { MDATA[c] = (int *) New_Track_Buffer(MTRACK[c]);
        if (n < MTRACK[c]->dmax)
          n = MTRACK[c]->dmax;
      }
    if (n > 0)
      MBUFFER = Malloc(sizeof(int64)*(n/sizeof(int)),"Allocating track buffer");
    else
      MBUFFER = NULL;
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

  schema = oneSchemaCreateFromText(One_Schema);
  file1  = oneFileOpenWriteNew("-",schema,"daz",true,1);
  oneAddProvenance(file1,Prog_Name,"1.0","%s >?.daz",command);

  //  Display each read (and/or QV streams) in the active DB according to the
  //    range pairs in pts[0..reps) and according to the display options.

  { DAZZ_READ  *reads;
    char       *read, *arrow, **entry;
    int         c, b, e, i, m;
    int         substr;
    int         map, last;
    char        qvname[5] = { 'D', 'C', 'I', 'M', 'S' };

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

    for (i = 0; i < MTOP; i++)
      { oneInt(file1,0) = i;
        oneWriteLine(file1,'X',strlen(MASK[i]),MASK[i]);
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

        if (map > 0 && findx[map-1] <= b && b < findx[map])
          ;
        else
          { map = 0;
            while (b >= findx[map])
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
            flags = r->flags;
            qv    = (flags & DB_QV);

            if (substr)
              { fst = iter->beg;
                lst = iter->end;
              }
            else
              { fst = 0;
                lst = len;
              }

            if (DOFLN && i >= findx[map])
              { if (strcmp(ffile[map+1],ffile[last]) != 0)
                  { oneInt(file1,0) = 0;
                    oneWriteLine(file1,'f',strlen(ffile[map+1]),ffile[map+1]);
                    last = map+1;
                  }
                if (!DOHDR || DAM)
                  map += 1;
              }

            Load_Read(db,i,read,1);

            oneInt(file1,0) = FirstRead + (i+1);
            oneWriteLine(file1,'R',lst-fst,read+fst);

            if (DOHDR)
              { if (DAM)
                  { char header[MAX_NAME];

                    FSEEKO(hdrs,r->coff,SEEK_SET)
                    FGETS(header,MAX_NAME,hdrs)
                    header[strlen(header)-1] = '\0';
                    oneWriteLine(file1,'H',strlen(header),header);
                  }
                else
                  { if (i >= findx[map])
                      { map += 1;
                        oneWriteLine(file1,'H',strlen(fhead[map]),fhead[map]);
                      }
                  }
              }

            if (DOWEL)
              { oneInt(file1,0) = r->origin;
                oneInt(file1,1) = r->fpulse;
                oneInt(file1,2) = r->fpulse+len;
                if (DAM)
                  oneWriteLine(file1,'G',0,NULL);
                else
                  oneWriteLine(file1,'W',0,NULL);
                if (Quiva_DB && qv > 0)
                  { oneInt(file1,0) = qv;
                    oneWriteLine(file1,'Q',0,NULL);
                  }
                else if (Arrow_DB)
                  { int   j, snr[4];
                    int64 big;

                    big = *((uint64 *) &(r->coff));
                    for (j = 0; j < 4; j++)
                      { snr[3-j] = (big & 0xffff);
                        big    >>= 16;
                      }
                    oneInt(file1,0) = snr[0];
                    oneInt(file1,1) = snr[1];
                    oneInt(file1,2) = snr[2];
                    oneInt(file1,3) = snr[3];
                    oneWriteLine(file1,'N',0,NULL);
                  }
              }

            for (m = 0; m < MTOP; m++)
              { int   *d, j;
                int64  f;

                d = MDATA[m];
                f = (Load_Track_Data(MTRACK[m],i,d) >> 2);

                for (j = 0; j < f; j++)
                  MBUFFER[j] = d[j];

                oneInt(file1,0) = m;
                oneWriteLine(file1,'T',f,MBUFFER);
              }

            if (DOARW)
              { Load_Arrow(db,i,arrow,1);
                oneWriteLine(file1,'A',lst-fst,arrow+fst);
              }

            if (DOQVS)
              { int k;

                Load_QVentry(db,i,entry,1);
                for (k = 0; k < 5; k++)
                  oneWriteLine(file1,qvname[k],lst-fst,entry[k]+fst);
              }
          }
      }
  }

  oneFileClose(file1);
  oneSchemaDestroy(schema);

  if (input_pts)
    { fclose(input);
      free(iter);
    }
  else
    free(pts);

  if (DAM)
    fclose(hdrs);
  else
    Free_DB_Stub(stub);
  Close_DB(db);

  free(command);

  exit (0);
}
