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
SpiceChar relate[SHORTLEN];
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
SpiceChar lmpoolString[LONGLEN];

LONGLONG nRows;
LONGLONG iRow;
LONGLONG iFacetNumber;
double dVNormal[3];

SpiceInt targetId;
SpiceInt nPool;
SpiceBoolean foundPool;

char* doDebug = getenv("GFTAG_DEBUG");

int iLoop;

int skipPoolFailPrint;

pGFTRN pGftRn = 0;

  do {

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

    /* Open normal vector FITS table */
    if (gftag_normals_open(normalsFITS,&pGftRn)) break;
    if (doDebug) { fprintf(stdout,"Opened normal vector FITS file [%s]\n",normalsFITS); }
    if (1>(nRows=gftag_normals_num_rows(pGftRn))) break;
    if (doDebug) { fprintf(stdout,"- Number of rows = %lld\n",nRows); }

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
    if (doDebug) { fprintf(stdout,"Loaded unit sphere radii [%s]\n", lmpoolString); }

    /* Loop over rows in the normal vectors table of the FITS file */
    for (iRow=0; iRow<nRows; ++iRow) {

      /* Initialize cnfine from cnfine0 */
      if (gftag_normals_next(pGftRn,&iFacetNumber,dVNormal)) {
        fprintf(stderr,"Failed to get facet number or normal for row [%lld]\n", iRow+1);
        break;
      }

      /* Get the surface point for this normal
       * - on a unit sphere, spoint is identical to unit normal vector at that point
       */
      vpack_c(dVNormal[0],dVNormal[1],dVNormal[2],spoint);
      vhat_c(spoint,spoint);

      /* Initialize cnfine from cnfine0; set relate to EQUALS */
      copy_c(&cnfine0,&cnfine);
      repmc_c("=","=","=",SHORTLEN,relate);


      /****************************************************************/
      /************* HERE'S THE BEEF **********************************/

      /* Find time(s) at which incidence angle equals that in refval */
      gfilum_c( method, angtype, target, illmn, fixref, abcorr, obsrvr, spoint
              , relate, refval, adjust, step, nintvls, &cnfine, &result);

      /****************************************************************/
      /****************************************************************/

      /* Handle error - go to next plate */
      if (failed_c()) {
        fprintf(stderr,"Failed glilum_c for row [%lld]\n", iRow+1);
        reset_c();
        continue;
      }

      /******************** TODO START ************************************/
      /* If cardinality of result is zero, look for local minimum instead */
      /******************** TODO END **************************************/

      /* Log cardinality of result, plus first and last times */
      if (doDebug) {
      long wncard = wncard_c(&result);
        fprintf(stdout,"gfilum_c found [%ld] solutions for facet [%lld]"
               , wncard, iFacetNumber);
        if (wncard>0) {
          timout_c(SPICE_CELL_ELEM_D(&result, 0), "YYYY-MM-DD/HR:MN:SC.### ::TDB", SHORTLEN, sTDB);
          fprintf(stdout,":  first TDB=%s",sTDB);
          if (wncard>1) {
            timout_c(SPICE_CELL_ELEM_D(&result, (wncard<<1)-1), "YYYY-MM-DD/HR:MN:SC.### ::TDB", SHORTLEN, sTDB);
            fprintf(stdout,"; last TDB=%s",sTDB);
          }
        }
        fprintf(stdout,".\n");

      } /* if (doDebug) */
    } /* for (iRow=...) */

    if (failed_c()) { break; }

  } while (0);

  kclear_c();
  if (pGftRn) { gftag_normals_close(&pGftRn); }
  
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

GFTAG_STEP_SECONDS       = 600.0
GFTAG_INCIDENCE_DEG      = 65.0

GFTAG_TDB_ENDPOINTS     += @2020-07-04-00:00:00
GFTAG_TDB_ENDPOINTS     += @2020-07-05-00:00:00


\begintext
#endif
