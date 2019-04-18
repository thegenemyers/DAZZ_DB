#include <stdlib.h>
#include <stdio.h>

#include "DB.h"

int main(int argc, char *argv[])
{ char   code, which;
  int64  total;
  int    len, rno, mno, qual;
  int    vec[4];
  char  *buffer[256];
  char **mname;
  int  **masks;

  (void) argv;

  if (argc > 1)
    { fprintf(stderr,"Usage: DBa2b <(ascii) >(binary)\n");
      exit (1);
    }

  while (scanf(" %c",&code) == 1)       //  Header lines
    if (code == '@' || code == '+')
      { scanf(" %c",&which);
        fwrite(&code,sizeof(char),1,stdout);
        fwrite(&which,sizeof(char),1,stdout);
        if (which == 'T')
          { scanf("%d %lld",&mno,&total);
            fwrite(&mno,sizeof(int),1,stdout);
            fwrite(&total,sizeof(int64),1,stdout);
            if (code == '@')
              { masks[mno] = (int *) malloc(sizeof(int)*2*total);
                scanf(" %d",&len);
                mname[mno] = (char *) malloc(sizeof(char)*(len+1));
                scanf(" %s",mname[mno]);
                fwrite(&len,sizeof(int),1,stdout);
                fwrite(mname[mno],sizeof(char),len,stdout);
              }
          }
        else
          { scanf(" %lld",&total);
            if (which == 'M')
              { masks = (int **) malloc(sizeof(int *)*total);
                mname = (char **) malloc(sizeof(char *)*total);
              }
            else if (code == '@')
              buffer[(int) which] = malloc(total+1);
            fwrite(&total,sizeof(int64),1,stdout);
          }
      }
    else
      { ungetc(code,stdin);
        break;
      }

  buffer['A'] = buffer['c'] = buffer['d'] = 
  buffer['i'] = buffer['m'] = buffer['s'] = buffer['S'];

  while (scanf(" %c",&code) == 1)       //  For each data line do
    { fwrite(&code,sizeof(char),1,stdout);
      switch (code)
      { case 'R':                         //  Read
          scanf(" %d",&rno);
          fwrite(&rno,sizeof(int),1,stdout);
          break;
        case 'Q':                         //  Read
          scanf(" %d",&qual);
          fwrite(&qual,sizeof(int),1,stdout);
          break;
        case 'L':                         //  Well, Pulse range
          scanf(" %d %d %d",vec,vec+1,vec+2);
          fwrite(vec,sizeof(int),3,stdout);
          break;
        case 'N':                         //  SNR values
          scanf(" %d %d %d %d",vec,vec+1,vec+2,vec+3);
          fwrite(vec,sizeof(int),4,stdout);
          break;
        case 'S': case 'A': case 'c':     //  DNA strings (2-bit compressible)
          scanf(" %d",&len);
          scanf(" %s",buffer[(int) code]);
          if (code == 'A')
            Number_Arrow(buffer[(int) code]);
          else
            Number_Read(buffer[(int) code]);
          Compress_Read(len,buffer[(int) code]);
          fwrite(&len,sizeof(int),1,stdout);
          fwrite(buffer[(int) code],sizeof(char),COMPRESSED_LEN(len),stdout);
          break;
        case 'H': case 'I': case 'F':     //  All other string fields
        case 'd': case 'i':
        case 'm': case 's':
          scanf(" %d",&len);
          scanf(" %s",buffer[(int) code]);
          fwrite(&len,sizeof(int),1,stdout);
          fwrite(buffer[(int) code],sizeof(char),len,stdout);
          break;
        case 'T':                         //  Mask
          scanf("%d %d",&mno,&len);
          for (int i = 0; i < len; i++)
            scanf(" %d %d",masks[mno]+2*i,masks[mno]+2*i+1);
          fwrite(&mno,sizeof(int),1,stdout);
          fwrite(&len,sizeof(int),1,stdout);
          fwrite(masks[mno],sizeof(int),2*len,stdout);
      }
    }

  exit (0);
}
