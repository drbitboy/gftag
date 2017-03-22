#include "gftag_normals.h"
#if 0
#define __GFTAG_NORMALS__
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"

typedef struct {
  fitsfile *fptr;
  long nRows;;
  long nextRow;
  int status;
  int iColXYZN[4];
  long iFacetNumber;
  double vNormal[3];
} GFTRN, *pGFTRN, **ppGFTRN;

int gftag_normals_open(char* filename,ppGFTRN ppGftRn);
int gftag_normals_close(ppGFTRN ppGftRn);
int gftag_normals_next(pGFTRN pGftRn);
int gftag_normals_set_row(pGFTRN pGftRn,long nextRow);
long gftag_normals_num_rows(pGFTRN pGftRn);

#endif /* __GFTAG_NORMALS__ */