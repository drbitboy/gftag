#include "gftag_normals.h"
#include "SpiceUsr.h"

int
main(int argc, char** argv) {
#define SHORTLEN 60
#define LONGLEN 2000
SpiceChar method[] = { "ELLIPSOID" };
SpiceChar angtype[] = { "INCIDENCE" };
SpiceChar target[SHORTLEN];
SpiceChar illmn[SHORTLEN];
SpiceChar fixref[SHORTLEN];
SpiceChar abcorr[SHORTLEN];
SpiceChar obsrvr[SHORTLEN];
SpiceDouble spoint[3];
SpiceChar relateLT[] = { "<" };
SpiceChar relateABSMIN[] = { "ABSMIN" };
SpiceChar relateABSMAX[] = { "ABSMAX" };
SpiceChar* relateUsed;
SpiceDouble refval;
SpiceDouble adjust = 0.;
SpiceDouble step;
#define MAXINTERVALS 100
#define MAXWIN (MAXINTERVALS<<1)
SpiceInt nintvls = MAXINTERVALS;
SPICEDOUBLE_CELL (cnfine0, MAXWIN);
SPICEDOUBLE_CELL (cnfine, MAXWIN);
SPICEDOUBLE_CELL (result, MAXWIN);
SpiceDouble et2[2];

SpiceChar sTDB[SHORTLEN];

SpiceChar k2l[LONGLEN];
SpiceChar normalsFITS[LONGLEN];
SpiceChar normalsFITSmd5[SHORTLEN];
SpiceChar commTiltFITSout[LONGLEN];
SpiceChar lmpoolString[LONGLEN];

long long nRows;
long long iRow;
long long debugRows = 0;
long long iFacetNumber;
double dVNormal[3];
double dUVSurfpt[3];

SpiceInt targetId;
SpiceInt nPool;
SpiceBoolean foundPool;

char* doDebug = getenv("GFTAG_DEBUG");

int iLoop;

int skipPoolFailPrint;

pGFTRN pGftRn = 0;

long wncard;

float* pEarthEmission = 0;
float* pSigma = 0;
SpiceDouble dpr = dpr_c();
SpiceDouble phaseAngle;
SpiceDouble incidenceAngle;
SpiceDouble emissionAngle;

fitsfile* inFptr = 0;
fitsfile* outFptr = 0;
int rtn;
int fitsStatus;
int* pFitsStatus = &fitsStatus;
int nCols;

  do {

    if (1!=sscanf(doDebug?doDebug:"","%Ld",&debugRows)) {
      debugRows = 0;
    }

    /* Ensure CSPICE errors do not abort this program */
    erract_c("set",0,"report");

    /* Command-line arguments are SPICE kernels */
    for (iLoop=1; !failed_c() && iLoop<argc; ++iLoop) {
      if (doDebug) fprintf(stdout,"FURNSHing [%s] ...",argv[iLoop]);
      furnsh_c(argv[iLoop]);
      if (doDebug) fprintf(stdout," ... done.\n");
    }
    if (failed_c()) break;

    /* Macro to read string Kernel Pool Variables (KPVs),
     * - With result checking
     * - The int skipPoolFailPrint determines if a failure is logged to stderr
     */
#   define GFT_GC(START,ROOM,SIN,SOUT) \
    gcpool_c(SIN,START,ROOM,sizeof(SOUT),&nPool,SOUT,&foundPool); \
    if (!(nPool==ROOM && foundPool) || failed_c()) { \
      if (!skipPoolFailPrint) { \
        fprintf(stderr,"Failed to get KPV [%s], START=%ld, ROOM=%ld\n" \
               , SIN, (long)START,(long)ROOM); \
      } \
      break; \
    } \
    if (doDebug) fprintf(stdout,"Got KPV [%s] = '%s'\n", SIN, SOUT)

    /* Get and FURNSH kernels in GFTAG_KERNELS_2_LOAD
     * - Loop until none are found
     * - Turn off logging
     ( - First missing item will break and exit this inner loop
     */
    skipPoolFailPrint = 1;

    /* Loop over this KPV until one is missing */
    for (iLoop=0; !failed_c() ; ++iLoop) {
      GFT_GC(iLoop,1,"GFTAG_KERNELS_2_LOAD",k2l);
      if (doDebug) fprintf(stdout,"FURNSHing [%s] ...",k2l);
      furnsh_c(k2l);
      if (doDebug) fprintf(stdout," ... done.\n");
    }
    if (failed_c()) break;

    /* Turn on logging */
    skipPoolFailPrint = 0;

    /* Get string KPVs; missing values will exit out do-while loop and end program */
    GFT_GC(0,1,"GFTAG_TARGET", target);
    GFT_GC(0,1,"GFTAG_TARGET_FIXREF", fixref);
    GFT_GC(0,1,"GFTAG_ILLUMINATOR", illmn);
    GFT_GC(0,1,"GFTAG_OBSERVER", obsrvr);
    GFT_GC(0,1,"GFTAG_ABCORR", abcorr);
    GFT_GC(0,1,"GFTAG_NORMALS_FITS", normalsFITS);
    GFT_GC(0,1,"GFTAG_NORMALS_FITS_MD5", normalsFITSmd5);
    GFT_GC(0,1,"GFTAG_COMM_TILT_FITS_OUT", commTiltFITSout);

    /* Open normal vector FITS table */
    if (gftag_normals_open(normalsFITS,&pGftRn)) break;
    if (doDebug) { fprintf(stdout,"Opened normal vector FITS file [%s]; MD5=[%s]\n",normalsFITS,normalsFITSmd5); }
    nRows = gftag_normals_num_rows(pGftRn);
    if (doDebug) { fprintf(stdout,"- Number of rows = %lld\n",nRows); }
    if (1>nRows) break;

    /* Allocate cleared memory for Earth emission angles and Sigma values */
    if (!(pEarthEmission=(float*)calloc(nRows, sizeof(float)))) {
      fprintf(stderr,"Failed to allocate memory for pEarthEmission\n");
      break;
    }

    if (!(pSigma=(float*)calloc(nRows, sizeof(float)))) {
      fprintf(stderr,"Failed to allocate memory for pSigma\n");
      break;
    }

    /* Macro for double-precision floating-point KPVs */
#   define GFT_GD(START,ROOM,SIN,DOUT) \
    gdpool_c(SIN,START,ROOM,&nPool,&DOUT,&foundPool); \
    if (!(nPool==ROOM && foundPool) || failed_c()) { \
      if (!skipPoolFailPrint) { \
        fprintf(stderr,"Failed to get KPV [%s], START=%ld, ROOM=%ld\n" \
               , SIN, (long)START,(long)ROOM); \
      } \
      break; \
    } \
    if (doDebug) fprintf(stdout,"Got KPV [%s] = %lg\n", SIN, DOUT)

    /* These two are required; convert incidence angle to radians */
    GFT_GD(0,1,"GFTAG_STEP_SECONDS",step);
    GFT_GD(0,1,"GFTAG_INCIDENCE_DEG",refval);
    refval *= rpd_c();

    /* Get TDB endpoint pairs, two at a time
     * - Loop until none are found
     * - Turn off logging
     * - First missing item will break and exit this inner loop
     */
    skipPoolFailPrint = 1;
    for (iLoop=0; !failed_c(); iLoop+=2) {
      /* Get two values */
      GFT_GD(iLoop,2,"GFTAG_TDB_ENDPOINTS",et2[0]);
      /* Append to cell cnfine0 */
      appndd_c(et2[0],&cnfine0);
      appndd_c(et2[1],&cnfine0);
      if (doDebug) fprintf(stdout,"Added TDB endpoints [%.3lf,%.3lf]; duration=%.3lf\n", et2[0], et2[1], et2[1]-et2[0]);
    }
    if (failed_c()) break;

    /* Validate pairs to make a window out of cnfine0 */
    wnvald_c(MAXWIN,iLoop,&cnfine0);
    if (failed_c()) break;

    /* Get target ID  from name */
    bods2c_c(target,&targetId,&foundPool);
    if (failed_c() || !foundPool) {
      fprintf(stderr,"Failed to get body ID for body name [%s]\n", target);
      break;
    }
    if (doDebug) { fprintf(stdout,"Target [%s] is body ID %ld\n", target, (long)targetId); }

    /* Build in-memory string to define a unit sphere shape model for the body,
     * to override any ellipsoid model that might have been already loaded via a PCK
     */
    repmi_c("BODY#_RADII = ( 1.0 1.0 1.0 )","#",targetId,LONGLEN,lmpoolString);
    lmpool_c(lmpoolString,LONGLEN,1);
    if (failed_c()) {
      fprintf(stderr,"Failed to get load radii via lmpool_c for body name/ID [%s/%ld]\n", target,(long)targetId);
      break;
    }
    if (doDebug) { fprintf(stdout,"Loaded unit sphere radii [%s]\n", lmpoolString); fflush(stdout); }

    /* Loop over rows in the normal vectors table of the FITS file */
    for (iRow=0; iRow<nRows && (!doDebug || debugRows<1 || iRow<debugRows); ++iRow) {

      pSigma[iRow] = pEarthEmission[iRow] = -1.0;

      /* Get next normal and surface point */
      if (gftag_normals_next(pGftRn,&iFacetNumber,dVNormal,0,dUVSurfpt)) {
        fprintf(stderr,"Failed to get facet number or normal for row [%lld]\n", iRow+1);
        break;
      }

      /* Covnert the surface point to a unit vector
       * - on a unit sphere, unit vector normal is identical to the surface point
       */
      vpack_c(dUVSurfpt[0],dUVSurfpt[1],dUVSurfpt[2],spoint);
      vhat_c(spoint,spoint);

      /* Initialize cnfine from cnfine0 */
      copy_c(&cnfine0,&cnfine);


      /****************************************************************/
      /************* HERE'S THE BEEF **********************************/

      /* Find time window(s) at which incidence angle is less than that
       * in refval; in nominal case, first and second values in first
       * window pair will be morning and afternoon, respectively, when
       * incidence angle is increasing and decreasing, respectively.
       * N.B. Edge cases
       *      1) No windows retuned in result:  incidence angle never
       *         crosses refval.  
       *      2) One window returned in result, and end of result window
       *         matches end of cnfine0 window, so incidence angle was
       *         less than refval at end of the returned window, and
       *         there is no afternoon event with incidence == refval.
       */
      gfilum_c( method, angtype, target, illmn, fixref, abcorr, obsrvr, spoint
              , relateLT, refval, adjust, step, nintvls, &cnfine, &result);

      /****************************************************************/
      /****************************************************************/

      /* Handle error - go to next plate */
      if (failed_c()) {
        fprintf(stderr,"Failed glilum_c for row [%lld]\n", iRow+1);
        reset_c();
        continue;
      }

      wncard = wncard_c(&result);
      relateUsed = relateLT;

      /* If cardinality of result is zero, look for absolute min instead */
      if (wncard == 0) {
        copy_c(&cnfine0,&cnfine);
        gfilum_c( method, angtype, target, illmn, fixref, abcorr, obsrvr, spoint
                , relateABSMIN, 0.0, 0.0, step, nintvls, &cnfine, &result);

        /* Handle error - go to next plate */
        if (failed_c()) {
          fprintf(stderr,"Failed ABSMIN glilum_c for row [%lld]\n", iRow+1);
          reset_c();
          continue;
        }

        wncard = wncard_c(&result);
        relateUsed = relateABSMIN;

        /* Get that incidence angle - assumes illumn is Sun */
        illum_c( target, SPICE_CELL_ELEM_D(&result,1), abcorr, obsrvr
               , spoint, &phaseAngle, &incidenceAngle, &emissionAngle);

        /* Handle error - go to next plate */
        if (failed_c()) {
          fprintf(stderr,"Failed illum_c for row [%lld]\n", iRow+1);
          reset_c();
          continue;
        }

        /* If ABSMIN incidence angle is less than the reference value,
         * then use the ABSMAX instead
         */
        if (incidenceAngle < refval) {
          copy_c(&cnfine0,&cnfine);
          gfilum_c( method, angtype, target, illmn, fixref, abcorr, obsrvr, spoint
                  , relateABSMAX, 0.0, 0.0, step, nintvls, &cnfine, &result);

          /* Handle error - go to next plate */
          if (failed_c()) {
            fprintf(stderr,"Failed ABSMAX gfilum_c for row [%lld]\n", iRow+1);
            reset_c();
            continue;
          }

          wncard = wncard_c(&result);
          relateUsed = relateABSMAX;

        }

      /* If cardinality of result is one, check for end of window */
      } else if (wncard == 1) {
        if ((et2[1]-SPICE_CELL_ELEM_D(&result,1)) < 1.0) {
          fprintf(stderr,"No afternoon event for row [%lld]\n", iRow+1);
          wncard = 0;
        }
      }

      if (wncard>0) {

        /* Get that incidence angle
         * - Assumes illumn is Sun
         */
        illum_c( target, SPICE_CELL_ELEM_D(&result,1), abcorr, obsrvr
               , spoint, &phaseAngle, &incidenceAngle, &emissionAngle);
        if (failed_c()) {
          reset_c();
          pSigma[iRow] = -1.0;
        } else {
          pSigma[iRow] = dpr * incidenceAngle;
        }

        /* Convert the normal to a unit vector
         * - on a unit sphere, spoint is identical to unit normal vector at that point
         */
        vpack_c(dVNormal[0],dVNormal[1],dVNormal[2],spoint);
        vhat_c(spoint,spoint);

        /* Get the emission angle wrt that normal */
        illum_c( target, SPICE_CELL_ELEM_D(&result,1), abcorr, obsrvr
               , spoint, &phaseAngle, &incidenceAngle, &emissionAngle);

        if (failed_c()) {
          pEarthEmission[iRow] = -1.0;
          reset_c();
        } else {
          pEarthEmission[iRow] = dpr * emissionAngle;
        }

      }

      /* Log cardinality of result, plus first and last times */
      if (doDebug) {
        fflush(stderr);
        fprintf(stdout,"gfilum_c found [%ld] [%-6s] solutions for facet [%lld]"
               , wncard, relateUsed, iFacetNumber);
        if (wncard>0) {
          timout_c(SPICE_CELL_ELEM_D(&result, 0), "YYYY-MM-DD/HR:MN:SC.### ::TDB", SHORTLEN, sTDB);
          fprintf(stdout,":  first TDB=%s",sTDB);
          timout_c(SPICE_CELL_ELEM_D(&result, 1), "YYYY-MM-DD/HR:MN:SC.### ::TDB", SHORTLEN, sTDB);
          fprintf(stdout,"; solution TDB=%s",sTDB);
        }
        fprintf(stdout,"; Earth emi=%.3lfdeg; Solar inc(sphere)=%.3lfdeg.\n",pEarthEmission[iRow],pSigma[iRow]);

        fflush(stdout);

      } /* if (doDebug) */
    } /* for (iRow=...) */

    if (failed_c()) { break; }

    if (gftag_normals_close(&pGftRn)) { break; }

    /* Open new FITS table */
    if ((rtn=fits_create_file(&outFptr,commTiltFITSout,pFitsStatus))) {
      fprintf(stderr,"--> Failed to open new FITS file [%s] in gftag\n\n",commTiltFITSout);
      fits_report_error(stderr,fitsStatus);
      break;
    }

    /* Open old normals FITS file */
    if ((rtn=fits_open_file(&inFptr,normalsFITS,READONLY,pFitsStatus))) {
      fprintf(stderr,"--> Failed to open normals FITS file [%s] in gftag\n\n",normalsFITS);
      fits_report_error(stderr,fitsStatus);
      break;
    }

    /* Copy normals FITS file to comm tilt FITS file */
    if ((rtn=fits_copy_file(inFptr,outFptr,1,1,1,pFitsStatus))) {
      fprintf(stderr,"--> Failed to copy normals FITS file [%s] to comm tilt map FITS file [%s] in gftag\n\n",normalsFITS,commTiltFITSout);
      fits_report_error(stderr,fitsStatus);
      break;
    }

    /* Move to PDU (first HDU) in comm tilt map FITS file */
    if ((rtn=fits_movabs_hdu(outFptr,1,0,pFitsStatus))) {
      fprintf(stderr,"--> Failed to move to PDU in normals FITS file [%s] in gftag\n\n",normalsFITS);
      fits_report_error(stderr,fitsStatus);
      break;
    }

    /* Add or update comm tilt map keywords to PDU */
    if (!rtn) { rtn=fits_update_key(outFptr,TSTRING,"MAP_NAME","COMM_TILT_MAP","",pFitsStatus); }
    if (!rtn) { rtn=fits_update_key(outFptr,TSTRING,"PRODNAME","comm_tilt_map.fits","",pFitsStatus); }
    timout_c(et2[0], "YYYY-MM-DD/HR:MN:SC.### TDB ::TDB", SHORTLEN, sTDB);
    if (failed_c()) { break; }
    if (!rtn) { rtn=fits_update_key(outFptr,TSTRING,"IN_DATE",sTDB,"",pFitsStatus); }
    if (!rtn) { rtn=fits_update_key(outFptr,TSTRING,"IN_VFILE",normalsFITS,"",pFitsStatus); }
    if (!rtn) { rtn=fits_update_key(outFptr,TSTRING,"IN_VMD5",normalsFITSmd5,"",pFitsStatus); }

    if (rtn) {
      fprintf(stderr,"--> Failed to update keywords in PDU for comm tilt map FITS file [%s] in gftag\n\n",commTiltFITSout);
      fits_report_error(stderr,fitsStatus);
      break;
    }

    /* Move to second HDU in comm tilt map  FITS file */
    if ((rtn=fits_movabs_hdu(outFptr,2,0,pFitsStatus))) {
      fprintf(stderr,"--> Failed to move to first HDU in normals FITS file [%s] in gftag\n\n",commTiltFITSout);
      fits_report_error(stderr,fitsStatus);
      break;
    }

    /* Get number of column in comm tilt map FITS file */
    if ((rtn=fits_get_num_cols(outFptr,&nCols,pFitsStatus))) {
      fprintf(stderr,"--> Failed to get number of columns in comm tilt map FITS file [%s] in gftag\n\n",commTiltFITSout);
      fits_report_error(stderr,fitsStatus);
      break;
    }

    /* Delete all columns beyond the fourth */
    while (nCols>4) {
      if ((rtn=fits_delete_col(outFptr,nCols,pFitsStatus))) {
        fprintf(stderr,"--> Failed to delete column [%d] in comm tilt map FITS file [%s] in gftag\n\n",nCols,commTiltFITSout);
        fits_report_error(stderr,fitsStatus);
        break;
      }
      --nCols;
    }

    if (nCols>4) { break; }

    /* Add columns for Earth emission angle and for SIGMA */
    {
    char tunitN[9];
    double zero = 0.0;
      rtn = fits_insert_col(outFptr,++nCols,"EARTH_EMISSION_ANGLE", "1E", pFitsStatus);
      sprintf(tunitN,"TUNIT%d",nCols);
      if (!rtn) { rtn=fits_update_key(outFptr,TSTRING,tunitN,"DEGREES","",pFitsStatus); }
      if (!rtn) { rtn=fits_write_col(outFptr,TFLOAT,nCols,1,1,nRows,pEarthEmission,pFitsStatus); }

      if (!rtn) { rtn=fits_insert_col(outFptr,++nCols,"SIGMA", "1E", pFitsStatus); }
      sprintf(tunitN,"TUNIT%d",nCols);
      if (!rtn) { rtn=fits_update_key(outFptr,TSTRING,tunitN,"DEGREES","",pFitsStatus); }
      sprintf(tunitN,"TZERO%d",nCols);
      if (!rtn) { rtn=fits_update_key(outFptr,TDOUBLE,tunitN,&zero,"",pFitsStatus); }
      sprintf(tunitN,"TSCAL%d",nCols);
      if (!rtn) { rtn=fits_update_key(outFptr,TDOUBLE,tunitN,&zero,"",pFitsStatus); }
      if (!rtn) { rtn=fits_write_col(outFptr,TFLOAT,nCols,1,1,nRows,pSigma,pFitsStatus); }

    }

    if (rtn) {
      fprintf(stderr,"--> Failed to add column [%d] in comm tilt map FITS file [%s] in gftag\n\n",nCols,commTiltFITSout);
      fits_report_error(stderr,fitsStatus);
      break;
    }

  } while (0);

  kclear_c();
  if (pEarthEmission) { free(pEarthEmission); }
  if (pSigma) { free(pSigma); }
  if (pGftRn) { gftag_normals_close(&pGftRn); }
  if (outFptr) { fits_close_file(outFptr,pFitsStatus); }
  if (inFptr) { fits_close_file(inFptr,pFitsStatus); }
  
  return 0;
}
#if 0
\begindata

GFTAG_KERNELS_2_LOAD    += 'kernels/mk.tm'
GFTAG_TARGET             = 'BENNU'
GFTAG_TARGET_FIXREF      = 'IAU_BENNU'
GFTAG_ILLUMINATOR        = 'SUN'
GFTAG_OBSERVER           = 'EARTH'
GFTAG_ABCORR             = 'LT'
GFTAG_NORMALS_FITS       = 'g_1254cm_tru_nvp_0000n00000_v100.fits'
GFTAG_NORMALS_FITS_MD5   = 'c316ef05ff2793e4df9004ad5804592e'

GFTAG_COMM_TILT_FITS_OUT = '!gftag.fits'

GFTAG_STEP_SECONDS       = 600.0
GFTAG_INCIDENCE_DEG      = 65.0

GFTAG_TDB_ENDPOINTS     += @2020-07-04-00:00:00
GFTAG_TDB_ENDPOINTS     += @2020-07-05-00:00:00


\begintext
#endif
