#include "gftag_normals.h"
#if 0
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
int gftag_normals_next(pGFTRN pGftRn);
int gftag_normals_set_row(pGFTRN pGftRn,LONGLONG nextRow);
LONGLONG gftag_normals_num_rows(pGFTRN pGftRn);

#endif /* __GFTAG_NORMALS__ */


/**********************************************************************/
/* Clean up allocated structure */
int
gftag_normals_cleanup(ppGFTRN ppGftRn) {
int rtn;
  if (!ppGftRn) {
    fprintf(stderr,"--> Null ppGFTRN pointer passed to gftag_normals_cleanup()\n\n");
    return -1;
  }
  if (!*ppGftRn) {
    fprintf(stderr,"--> Null pGFTRN pointer passed to gftag_normals_cleanup()\n\n");
    return -2;
  }
  free(*ppGftRn);
  *ppGftRn = 0;
  return 0;
}


/**********************************************************************/
/* Close FITS file, clean up */
int
gftag_normals_close(ppGFTRN ppGftRn) {
int rtn;
  if (!ppGftRn) {
    fprintf(stderr,"--> Null ppGFTRN pointer passed to gftag_normals_close()\n\n");
    return -1;
  }
  if (!*ppGftRn) {
    fprintf(stderr,"--> Null pGFTRN pointer passed to gftag_normals_close()\n\n");
    return -2;
  }
  if ((*ppGftRn)->fptr) {
    fprintf(stderr,"--> Null (fitsfile*) pointer passed to gftag_normals_close()\n\n");
    return -3;
  }

  if ((rtn=fits_close_file((*ppGftRn)->fptr,&(*ppGftRn)->status))) {
    fprintf(stderr,"--> Failed to close FITS file in gftag_normals_close()\n\n");
    return -4;
  }
  return gftag_normals_cleanup(ppGftRn);
}


/**********************************************************************/
/* Close FITS file, clean up */
int gftag_normals_set_row(pGFTRN pGftRn,LONGLONG nextRow) {
  if (!pGftRn) {
    fprintf(stderr,"--> Null pGFTRN pointer passed to gftag_normals_set_row()\n\n");
    return -1;
  }
  pGftRn->nextRow = nextRow;
  return 0;
}


/**********************************************************************/
/* Allocat GFTRN structure, open FITS file to first table, confirm it
 * is normal vector table
 */
int gftag_normals_open(char* filename,ppGFTRN ppGftRn) {
pGFTRN pGftRn;
int rtn;
int myRtn = 0;
int iXYZN;
  /* Check input pointers */
  if (!filename) {
    fprintf(stderr,"\n--> Null filename passed to gftag_normals_open()\n\n");
    return -1;
  }
  if (!ppGftRn) {
    fprintf(stderr,"\n--> Null ppGFTRN pointer passed to gftag_normals_open() when trying to open FITS file [%s]\n\n",filename);
    return -2;
  }

  /* Allocate GFRTN struct filled with zeros */
  if (!(pGftRn = *ppGftRn = (pGFTRN) calloc(1,sizeof(GFTRN)))) {
    fprintf(stderr,"\n--> Failed to calloc memory for GFTRN struct in gftag_normals_open() when trying to open FITS file [%s]\n\n",filename);
    return -3;
  }

  /* Open FITS file at first table */
  if ((rtn=fits_open_table(&pGftRn->fptr,filename,READONLY,&pGftRn->status))) {
    fits_report_error(stderr,pGftRn->status);
    fprintf(stderr,"\n--> Failed to open FITS file [%s] in gftag_normals_open()\n\n",filename);
    gftag_normals_cleanup(ppGftRn);
    return -4;
  }

  /* Get number of rows */
  if ((rtn=fits_get_num_rowsll(pGftRn->fptr,&pGftRn->nRows,&pGftRn->status))) {
    fits_report_error(stderr,pGftRn->status);
    fprintf(stderr,"\n--> Failed to open FITS file [%s] in gftag_normals_open()\n\n",filename);
    return -5;
  }

  /* Get column number for columns named NORMAL VECTOR_X/Y/Z and FACET NUMBER */
  for (iXYZN=0; !rtn && iXYZN<4; ++iXYZN) {
  static char* colnames[4] = { "NORMAL VECTOR_X", "NORMAL VECTOR_Y", "NORMAL VECTOR_Z", "FACET NUMBER" };
    if ((rtn=fits_get_colnum(pGftRn->fptr,CASESEN,colnames[iXYZN],pGftRn->iColXYZN+iXYZN,&pGftRn->status))) {
      fits_report_error(stderr,pGftRn->status);
      fprintf(stderr,"\n--> Failed to get column number for [%s] in FITS file [%s] in gftag_normals_open()\n\n",colnames[iXYZN],filename);
      myRtn = -6 - iXYZN;
    }
  }

  if (!rtn) {
    if ((rtn=gftag_normals_set_row(*ppGftRn,1LL))) {
      fprintf(stderr,"\n--> Failed to set row number to %lld in FITS file [%s] in gftag_normals_open()\n\n",1LL,filename);
      myRtn = -10;
    }
  }

  if (rtn) {
    if ((rtn=gftag_normals_close(ppGftRn))) {
      fprintf(stderr,"\n--> Failed to close FITS file [%s] in gftag_normals_open()\n\n",filename);
      myRtn = -10;
    }
  }

  return myRtn;
}


/**********************************************************************/
/* Get normal vector components and facet number from next row */
int gftag_normals_next(pGFTRN pGftRn) {
  return 0;
}


/**********************************************************************/
/* Close FITS file, clean up */
LONGLONG gftag_normals_num_rows(pGFTRN pGftRn) {
  return 0;
}
