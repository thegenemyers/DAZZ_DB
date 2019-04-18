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
    { fprintf(stderr,"Usage: DBb2a <(binary) >(ascii)\n");
      exit (1);
    }

  if (fread(&code,sizeof(char),1,stdin) == 0)
    code = 0;
    
  while (code == '@' || code == '+')
    { fread(&which,sizeof(char),1,stdin);
      printf("%c %c",code,which);
      if (which == 'T')
        { fread(&mno,sizeof(int),1,stdin);
          fread(&total,sizeof(int64),1,stdin);
          printf("%d %lld",mno,total);
          if (code == '@')
            { masks[mno] = (int *) malloc(sizeof(int)*2*total);
              fread(&len,sizeof(int),1,stdin);
              mname[mno] = (char *) malloc(sizeof(char)*(len+1));
              fread(mname[mno],sizeof(char),len,stdin);
              printf(" %d %.*s",len,len,mname[mno]);
            }
          printf("\n");
        }
      else
        { fread(&total,sizeof(int64),1,stdin);
          if (which == 'M')
            { masks = (int **) malloc(sizeof(int *)*total);
              mname = (char **) malloc(sizeof(char *)*total);
            }
          else if (code == '@')
            buffer[(int) which] = malloc(total+1);
          printf(" %lld\n",total);
        }

      if (fread(&code,sizeof(char),1,stdin) == 0)
       code = 0;
    }

  buffer['A'] = buffer['c'] = buffer['d'] = 
  buffer['i'] = buffer['m'] = buffer['s'] = buffer['S'];

  while (code != 0)       //  For each data line do
    { switch (code)
      { case 'R':                         //  Read
          fread(&rno,sizeof(int),1,stdin);
          printf("R %d\n",rno);
          break;
        case 'Q':                         //  Qual Value
          fread(&qual,sizeof(int),1,stdin);
          printf("Q %d\n",qual);
          break;
        case 'L':                         //  Well, Pulse range
          fread(vec,sizeof(int),3,stdin);
          printf("L %d %d %d\n",vec[0],vec[1],vec[2]);
          break;
        case 'N':                         //  SNR values
          fread(vec,sizeof(int),4,stdin);
          printf("N %d %d %d %d\n",vec[0],vec[1],vec[2],vec[3]);
          break;
        case 'S': case 'A': case 'c':     //  DNA strings (2-bit compressible)
          fread(&len,sizeof(int),1,stdin);
          fread(buffer[(int) code],sizeof(char),COMPRESSED_LEN(len),stdin);
          Uncompress_Read(len,buffer[(int) code]);
          if (code == 'A')
            Letter_Arrow(buffer[(int) code]);
          else
            Lower_Read(buffer[(int) code]);
          printf("%c %d %.*s\n",code,len,len,buffer[(int) code]);
          break;
        case 'H': case 'I': case 'F':     //  All other string fields
        case 'd': case 'i':
        case 'm': case 's':
          fread(&len,sizeof(int),1,stdin);
          fread(buffer[(int) code],sizeof(char),len,stdin);
          printf("%c %d %.*s\n",code,len,len,buffer[(int) code]);
          break;
        case 'T':                         //  Mask
          fread(&mno,sizeof(int),1,stdin);
          fread(&len,sizeof(int),1,stdin);
          fread(masks[mno],sizeof(int),2*len,stdin);
          printf("T%d %d",mno,len);
          for (int i = 0; i < len; i++)
            printf(" %d %d",masks[mno][2*i],masks[mno][2*i+1]);
          printf("\n");
      }

      if (fread(&code,sizeof(char),1,stdin) == 0)
       code = 0;
    }

  exit (0);
}
