/********************************************************************************************
 *
 *  Concate in block order all "block tracks" <DB>.<track>.# into a single track
 *    <DB>.<track>
 *
 *  Author:  Gene Myers
 *  Date  :  June 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-vfd] <path:db|dam> <track:name> ...";

int main(int argc, char *argv[])
{ char *prefix;
  int   nblocks;
  FILE *aout, *dout;
  int   c;

  int   VERBOSE;
  int   FORCE;
  int   DELETE;

  //  Process arguments

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("Catrack")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vfd") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    FORCE   = flags['f'];
    DELETE  = flags['d'];

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"   -v: verbose\n");
        fprintf(stderr,"   -d: delete individual blocks after a successful concatenation\n");
        fprintf(stderr,"   -f: force overwrite of track if already present\n");
        exit (1);
      }
  }

  //  Open DB stub and get number of blocks

  { char *pwd, *root;
    int   i, plen, index, isdam;
    FILE *dstub;
    char *dstub_name;

    plen = strlen(argv[1]);
    if (strcmp(argv[1]+(plen-3),".dam") == 0)
      root = Root(argv[1],".dam");
    else
      root = Root(argv[1],".db");
    pwd = PathTo(argv[1]);
    prefix = Strdup(Catenate(pwd,PATHSEP,root,"."),"Allocating track name");

    dstub = fopen(Catenate(pwd,"/",root,".db"),"r");
    isdam = 0;
    if (dstub == NULL)
      { dstub = fopen(Catenate(pwd,"/",root,".dam"),"r");
        isdam = 1;
        if (dstub == NULL)
          { fprintf(stderr,"%s: Cannot find %s either as a .db or a .dam\n",Prog_Name,root);
            exit (1);
          }
      }
    dstub_name = Strdup(Catenate(pwd,"/",root,isdam?".dam":".db"),"Allocating db file name");
    if (dstub_name == NULL)
      exit (1);
    
    FSCANF(dstub,DB_NFILE,&nblocks)
    
    for (i = 0; i < nblocks; i++)
      { char prolog[MAX_NAME], fname[MAX_NAME];
        
        FSCANF(dstub,DB_FDATA,&index,fname,prolog)
      }
    
    FSCANF(dstub,DB_NBLOCK,&nblocks)

    fclose(dstub);
    free(dstub_name);
    free(pwd);
    free(root);
  }

  //  For each track do

  for (c = 2; c < argc; c++)
    { int   nfiles;
      int   tracktot, tracksiz;
      int64 trackoff;
      char  data[1024];
      void *anno;
      FILE *lfile = NULL;
      DAZZ_EXTRA *extra;
      int         nextra;
      int64       extail;

      //  Open the output .anno for writing, and output header stub

      if (VERBOSE)
        { fprintf(stderr,"\nConstructing %s%s:\n",prefix,argv[c]);
          fflush(stderr);
        }

      aout = fopen(Catenate(prefix,argv[c],".","anno"),"r");
      if (aout != NULL)
        { if (!FORCE)
            { fprintf(stderr,"%s: Track file %s%s.anno already exists!\n",Prog_Name,prefix,argv[c]);
              exit (1);
            }
          fclose(aout);
        }

      dout = fopen(Catenate(prefix,argv[c],".","data"),"r");
      if (dout != NULL)
        { if (!FORCE)
            { fprintf(stderr,"%s: Track file %s%s.data already exists!\n",Prog_Name,prefix,argv[c]);
              exit (1);
            }
          fclose(dout);
        }

      aout = Fopen(Catenate(prefix,argv[c],".","anno"),"w");
      if (aout == NULL)
        exit (1);
      dout = NULL;

      extra    = NULL;
      anno     = NULL;
      trackoff = 0;
      tracktot = tracksiz = 0;
      if (fwrite(&tracktot,sizeof(int),1,aout) != 1)
        SYSTEM_WRITE_ERROR
      if (fwrite(&tracksiz,sizeof(int),1,aout) != 1)
        SYSTEM_WRITE_ERROR

      //  OPen and catenate in each block .anno and .data file

      nextra = 0;
      nfiles = 0;
      while (1)
        { FILE *dfile, *afile;
          char *dfile_name, *afile_name;
          int   i, apos, size, esize, tracklen;
  
          afile_name = Strdup(Numbered_Suffix(prefix,nfiles+1,Catenate(".",argv[c],".","anno")),
                              "Allocating .anno file name");
          dfile_name = Strdup(Numbered_Suffix(prefix,nfiles+1,Catenate(".",argv[c],".","data")),
                              "Allocating .data file name");
          if (afile_name == NULL || dfile_name == NULL)
            goto error;
            
  
          afile = fopen(afile_name,"r");
          if (afile == NULL)
            break;
          dfile = fopen(Numbered_Suffix(prefix,nfiles+1,Catenate(".",argv[c],".","data")),"r");
          if (dfile == NULL && errno != ENOENT)
            { fprintf(stderr,"%s: The file %s is corrupted\n",Prog_Name,dfile_name);
              goto error;
            }
  
          if (nfiles > 0)
            fclose(lfile);
          lfile = afile;
  
          if (VERBOSE)
            { fprintf(stderr,"  Concatenating %s%d.%s ...\n",prefix,nfiles+1,argv[c]);
              fflush(stderr);
            }
    
          FFREAD(&tracklen,sizeof(int),1,afile)
          FFREAD(&size,sizeof(int),1,afile)
          if (size == 0)
            esize = 8;
          else
            esize = size;
  
          if (nfiles == 0)
            { tracksiz = size;
              if (dfile != NULL)
                { dout = Fopen(Catenate(prefix,argv[c],".","data"),"w");
                  if (dout == NULL)
                    goto error;
                }
              else
                { anno = Malloc(esize,"Allocating annotation record");
                  if (anno == NULL)
                    goto error;
                }
            }
          else
            { int escape = 1;
              if (tracksiz != size)
                { fprintf(stderr,"%s: Track block %d does not have the same annotation size (%d)",
                                 Prog_Name,nfiles+1,size);
                  fprintf(stderr," as previous blocks (%d)\n",tracksiz);
                }
              else if (dfile == NULL && dout != NULL)
                fprintf(stderr,"%s: Track block %d does not have data but previous blocks do\n",
                               Prog_Name,nfiles+1);
              else if (dfile != NULL && dout == NULL)
                fprintf(stderr,"%s: Track block %d has data but previous blocks do not\n",
                               Prog_Name,nfiles+1);
              else
                 escape = 0;
              if (escape)
                goto error;
            }
    
          if (dfile != NULL)
            { int64 dlen, d;
  
              if (esize == 4)
                { int anno4;
    
                  for (i = 0; i < tracklen; i++)
                    { FFREAD(&anno4,sizeof(int),1,afile)
                      anno4 += trackoff;
                      FFWRITE(&anno4,sizeof(int),1,aout)
                    }
                  FFREAD(&anno4,sizeof(int),1,afile)
                  dlen = anno4;
                }
              else
                { int64 anno8;
    
                  for (i = 0; i < tracklen; i++)
                    { FFREAD(&anno8,sizeof(int64),1,afile)
                      anno8 += trackoff;
                      FFWRITE(&anno8,sizeof(int64),1,aout)
                    }
                  FFREAD(&anno8,sizeof(int64),1,afile)
                  dlen = anno8;
                }
              trackoff += dlen;

              for (d = 1024; d < dlen; d += 1024)
                { FFREAD(data,1024,1,dfile)
                  FFWRITE(data,1024,1,dout)
                }
              d -= 1024;
              if (d < dlen)
                { FFREAD(data,dlen-d,1,dfile)
                  FFWRITE(data,dlen-d,1,dout)
                }
            }
          else
            { for (i = 0; i < tracklen; i++)
                { FFREAD(anno,esize,1,afile)
                  FFWRITE(anno,esize,1,aout)
                }
            }
  
          FSEEKO(afile,0,SEEK_END)
          FTELLO(apos,afile)
          if (dfile != NULL)
            extail = apos - (esize*(tracklen+1) + 2*sizeof(int)); 
          else
            extail = apos - (esize*tracklen + 2*sizeof(int)); 
          FSEEKO(afile,-extail,SEEK_END)

          if (extail >= 20)
            { if (nfiles == 0)
                { nextra = 0;
                  while (1)
                    if (Read_Extra(afile,afile_name,NULL))
                      break;
                    else
                      nextra += 1;
  
                  extra = (DAZZ_EXTRA *) Malloc(sizeof(DAZZ_EXTRA)*(nextra+1),"Allocating extras");
                  if (extra == NULL)
                    goto error;
                  FSEEKO(afile,-extail,SEEK_END)
  
                  for (i = 0; i < nextra; i++)
                    { extra[i].nelem = 0;
                      Read_Extra(afile,afile_name,extra+i);
                    }
                }
  
              else
                { for (i = 0; i < nextra; i++)
                    if (Read_Extra(afile,afile_name,extra+i))
                      { fprintf(stderr,"%s: File %s has fewer extras than previous .anno files\n",
                                       Prog_Name,afile_name);
                        goto error;
                      }
                  if (Read_Extra(afile,afile_name,extra+nextra) == 0)
                    { fprintf(stderr,"%s: File %s has more extras than previous .anno files\n",
                                     Prog_Name,afile_name);
                      goto error;
                    }
                }
            }
    
          tracktot += tracklen;
          nfiles   += 1;
          if (dfile != NULL)
            fclose(dfile);
        }

      if (nfiles == 0)
        { fprintf(stderr,"%s: Couldn't find first track block %s1.%s.anno\n",
                         Prog_Name,prefix,argv[c]);
          goto error;
        }
      else
        { char *byte;

          if (dout != NULL)
            { if (tracksiz == 4)
                { int anno4 = trackoff;
                  FFWRITE(&anno4,sizeof(int),1,aout)
                }
              else
                { int64 anno8 = trackoff;
                  FFWRITE(&anno8,sizeof(int64),1,aout)
                }
            }
  
          if (nextra == 0)
            { while (fread(&byte,1,1,lfile) == 1)
                FFWRITE(&byte,1,1,aout)
            }
          else
            { int i;
              for (i = 0; i < nextra; i++)
                Write_Extra(aout,extra+i);
            }
          fclose(lfile);

          FSEEKO(aout,0,SEEK_SET)
          FFWRITE(&tracktot,sizeof(int),1,aout)
          FFWRITE(&tracksiz,sizeof(int),1,aout)
        }

      if (nfiles != nblocks)
        { fprintf(stderr,"%s: Did not catenate all tracks of DB (nfiles %d != nblocks %d)\n",
	                 Prog_Name, nfiles, nblocks);
          goto error;
        }
  
      FCLOSE(aout);
      if (dout != NULL)
        FCLOSE(dout);

      if (DELETE)
        { int   i;
          char *name;

          for (i = 1; i <= nblocks ;i++)
            { name = Numbered_Suffix(prefix,i,Catenate(".",argv[c],".","anno"));
              if (unlink(name) != 0)
                fprintf(stderr,"%s: [WARNING] Couldn't delete file %s\n",Prog_Name,name);
              if (dout != NULL)
                { name = Numbered_Suffix(prefix,i,Catenate(".",argv[c],".","data"));
                  if (unlink(name) != 0)
                    fprintf(stderr,"%s: [WARNING] Couldn't delete file %s\n",Prog_Name,name);
                }
            }
        }
    
    }

  free(prefix);
  exit (0);

error:
  { char *name;

    fclose(aout);
    name = Catenate(prefix,argv[c],".","anno");
    if (unlink(name) != 0)
      fprintf(stderr,"%s: [WARNING] Couldn't delete file %s during abort\n",Prog_Name,name);
    if (dout != NULL)
      { fclose(dout);
        name = Catenate(prefix,argv[c],".","data");
        if (unlink(name) != 0)
          fprintf(stderr,"%s: [WARNING] Couldn't delete file %s during abort\n",Prog_Name,name);
      }
  }

  free(prefix);
  exit (1);
}
