#include <math.h>
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
  int iColLLRXYZN[7];
  long iFacetNumber;
  double vLatLonR[3];
  double vNormal[3];
} GFTRN, *pGFTRN, **ppGFTRN;

int gftag_normals_close(ppGFTRN ppGftRn);
int gftag_normals_open(char* filename,ppGFTRN ppGftRn);
int gftag_normals_next(pGFTRN pGftRn,LONGLONG* pIFacetNumber, double* vNormal, double* pLatLonR, double* pUvSurfpt);
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
  if (!(*ppGftRn)->fptr) {
    fprintf(stderr,"--> Null (fitsfile*) pointer passed to gftag_normals_close()\n\n");
    return -3;
  }

  if ((rtn=fits_close_file((*ppGftRn)->fptr,&(*ppGftRn)->status))) {
    fprintf(stderr,"--> Failed to close FITS file in gftag_normals_close()\n\n");
    fits_report_error(stderr,(*ppGftRn)->status);
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
int iLLRXYZN;
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
  for (iLLRXYZN=0; !rtn && iLLRXYZN<7; ++iLLRXYZN) {
  static char* colnames[7] = { "NORMAL VECTOR_X", "NORMAL VECTOR_Y", "NORMAL VECTOR_Z", "LATITUDE", "LONGITUDE", "RADIUS", "FACET NUMBER" };
    if ((rtn=fits_get_colnum(pGftRn->fptr,CASESEN,colnames[iLLRXYZN],pGftRn->iColLLRXYZN+iLLRXYZN,&pGftRn->status))) {
      fits_report_error(stderr,pGftRn->status);
      fprintf(stderr,"\n--> Failed to get column number for [%s] in FITS file [%s] in gftag_normals_open()\n\n",colnames[iLLRXYZN],filename);
      myRtn = -6 - iLLRXYZN;
    }
  }

  if (!rtn) {
    if ((rtn=gftag_normals_set_row(*ppGftRn,1LL))) {
      fprintf(stderr,"\n--> Failed to set row number to %lld in FITS file [%s] in gftag_normals_open()\n\n",1LL,filename);
      myRtn = -13;
    }
  }

  if (rtn) {
    if ((rtn=gftag_normals_close(ppGftRn))) {
      fprintf(stderr,"\n--> Failed to close FITS file [%s] in gftag_normals_open()\n\n",filename);
      myRtn = -14;
    }
  }

  return myRtn;
}


/**********************************************************************/
/* Close FITS file, clean up */
LONGLONG gftag_normals_num_rows(pGFTRN pGftRn) {
  if (!pGftRn) {
    fprintf(stderr,"--> Null pGFTRN pointer passed to gftag_normals_num_rows()\n\n");
    return -1LL;
  }
  return pGftRn->nRows;
}


/**********************************************************************/
/* Get normal vector components and facet number from next row */
int
gftag_normals_next(pGFTRN pGftRn,LONGLONG* pIFacetNumber, double* pVNormal, double* pLatLonR, double* pUvSurfpt) {
int rtn = 0;
int myRtn = 0;
const double nulval = -999.0;
int anynul;

  do {
  int iLLRXYZN;
  LONGLONG saveRow;

    if (!pGftRn) {
      fprintf(stderr,"--> Null pGFTRN pointer passed to gftag_normals_next()\n\n");
      myRtn = -1;
      break;
    }
    if (!pGftRn->fptr) {
      fprintf(stderr,"--> Null fitfile* pointer passed to gftag_normals_next()\n\n");
      myRtn = -2;
      break;
    }

    saveRow = pGftRn->nextRow;

    for (iLLRXYZN=0; !myRtn && !rtn && iLLRXYZN<7; ++iLLRXYZN) {
      if ((rtn=fits_read_col(pGftRn->fptr                     /* filelist* fptr       */
               , iLLRXYZN<6 ? TDOUBLE : TLONGLONG             /* datatype             */
               , pGftRn->iColLLRXYZN[iLLRXYZN]                /* colnum               */
               , saveRow                                      /* firstrow             */
               , 1LL, 1LL                                     /* firstelem, nelements */
               , (void*) &nulval                              /* DTYPE* nulval        */
               , iLLRXYZN<3
                 ? (void*) (pGftRn->vNormal+iLLRXYZN)         /* DTYPE* array         */
                 : iLLRXYZN<6                                 /* DTYPE* array         */
                   ? (void*) (pGftRn->vLatLonR+(iLLRXYZN-3))  /* DTYPE* array         */
                   : (void*) (&pGftRn->iFacetNumber)          /* DTYPE* array         */
               , &anynul                                      /* anynul               */
               , &pGftRn->status                              /* firstrow             */
               ))) {
        fits_report_error(stderr,pGftRn->status);
        fprintf(stderr,"--> fits_read_col() failed for column [%d] of row [%lld] in gftag_normals_next()\n\n"
                      , pGftRn->iColLLRXYZN[iLLRXYZN], saveRow);
        myRtn = -3;
      }
    }

    if (!myRtn && !rtn) {
      /* Copy successful results to non-null arguments */
      if (pIFacetNumber) { *pIFacetNumber = pGftRn->iFacetNumber; }
      if (pVNormal) { memcpy( pVNormal, pGftRn->vNormal, 3 * sizeof(double)); }
      if (pLatLonR) { memcpy( pLatLonR, pGftRn->vLatLonR, 3 * sizeof(double)); }
      if (pUvSurfpt) {
      static double rpd = M_PI / 180.0;
      double coslat = cos(pGftRn->vLatLonR[0] * rpd);
      double sinlat = sin(pGftRn->vLatLonR[0] * rpd);
      double coslon = cos(pGftRn->vLatLonR[1] * rpd);
      double sinlon = sin(pGftRn->vLatLonR[1] * rpd);
        pUvSurfpt[0] = coslat * coslon;
        pUvSurfpt[1] = coslat * sinlon;
        pUvSurfpt[2] = sinlat;
      }

      /* Increment row */
      ++pGftRn->nextRow;
    }
  } while(0);

  return myRtn;
}

/**********************************************************************/
/**********************************************************************/
/* Test code */
int
gftag_normals_test(int argc,char** argv) {
pGFTRN pGftRn = 0;
int rtn = 0;
  do {
  LONGLONG nRows;
  LONGLONG thisRow;
  LONGLONG thisFacetNumber;
  double normUVector[3];
  double llrVector[3];
  double surfptUVector[3];
    if (argc<2) {
      fprintf(stderr, "Usage:  %s normals_table.fit\n\n", argc>0 ? argv[0] : "test_gftag_normals");
      rtn = -1;
      break;
    }
    if ((rtn=gftag_normals_open(argv[1],&pGftRn))) {
      fprintf(stderr, "\n--> gftag_normals_open() failed for FITS filename [%s]\n\n",argv[1]);
      rtn = -2;
      break;
    }

    if (0>=(nRows=gftag_normals_num_rows(pGftRn))) {
      fprintf(stderr, "\n--> gftag_normals_num_rows() returned a non-positive value for FITS filename [%s]\n\n",argv[1]);
      rtn = -3;
      break;
    }
    fprintf(stdout, "There are %lld rows in the normal vectors table in FITS file [%s]\n",nRows,argv[1]);

    for (thisRow=1; !rtn && thisRow<=nRows; thisRow+=(nRows-1)) {
    LONGLONG saveRow;
      saveRow = pGftRn->nextRow;
      if (saveRow!=thisRow) {
        fprintf(stderr, "\n--> gftag_normals_test() failed with expected row mismatch (%lld!=%lld in FITS filename [%s]\n\n",saveRow,thisRow,argv[1]);
        rtn = -4;
      } else  if ((rtn=gftag_normals_next(pGftRn,&thisFacetNumber,normUVector,llrVector,surfptUVector))) {
        fprintf(stderr, "\n--> gftag_normals_next() failed for row [%lld] in FITS filename [%s]\n\n",pGftRn->nextRow,argv[1]);
        rtn = -5 - (thisRow==1 ? 0 : 1) ;
      } else {
        fprintf(stdout, "\n%s row of %s:\n\
        facet number =  %10lld\n\
       normal vector = [%10.6lf %10.6lf %10.6lf]\n\
    Lat,ELong,Radius = [%10.5lf %10.5lf %10.5lf]\n\
  surfpt unit vector = [%10.6lf %10.6lf %10.6lf]\n\
"
                      , thisRow==1 ? "First" : "Last"
                      , argv[1]
                      , thisFacetNumber
                      , normUVector[0], normUVector[1], normUVector[2]
                      , llrVector[0], llrVector[1], llrVector[2]
                      , surfptUVector[0], surfptUVector[1], surfptUVector[2]
                      );
        if ((rtn=gftag_normals_set_row(pGftRn,saveRow + (nRows-1)))) {
          fprintf(stderr, "\n--> gftag_normals_set_row() failed for row [%lld] in FITS filename [%s]\n\n",saveRow+(nRows-1),argv[1]);
          rtn = -7;
        }
      } /* if ... else if ... else */
    } /* for (thisRow...) */
  } while (0);

  if (pGftRn) {
  int rtn2;
    if ((rtn2=gftag_normals_close(&pGftRn))) {
      fprintf(stderr, "\n--> gftag_normals_close() failed for FITS filename [%s]\n\n",argv[1]);
      rtn = -6;
    }
  }

  return rtn;
}

#ifdef __DO_MAIN_TEST__
int
main(int argc, char** argv) {
  return gftag_normals_test(argc,argv);
}
#endif
