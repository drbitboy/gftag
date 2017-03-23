#ifndef __GFTAG_NORMALS__
#define __GFTAG_NORMALS__
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"

typedef struct {
  fitsfile *fptr;
  LONGLONG nRows;;
  LONGLONG nextRow;
  int status;
  int iColXYZN[4];
  long iFacetNumber;
  double vNormal[3];
} GFTRN, *pGFTRN, **ppGFTRN;

int gftag_normals_close(ppGFTRN ppGftRn);
int gftag_normals_open(char* filename,ppGFTRN ppGftRn);
int gftag_normals_next(pGFTRN pGftRn,LONGLONG* pIFacetNumber, double* vNormal);
int gftag_normals_set_row(pGFTRN pGftRn,LONGLONG nextRow);
LONGLONG gftag_normals_num_rows(pGFTRN pGftRn);

#endif /* __GFTAG_NORMALS__ */
