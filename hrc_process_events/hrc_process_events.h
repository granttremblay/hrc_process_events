/*                                                                
**  Copyright (C) 1996-2010,2012  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */


/************************************************************************
 *JCC(5/1/00)-updated for tap corrections to include 'tap_ring_defs.h' and 
 *        to add tapfile, amps_3RD_raw and amps_tap_flag for EVENT_REC_P_T.
 *JCC(5/11/00)-compute corrections using the 'double' data type for the
 *             amplitudes; add amps_dd and amps_sh ; remove amps_3RD_raw ;
 *JCC(5/12/00)-saturation coefficients are a function of amp_sf;
 *             update open_amp_saturation_file and check_amp_saturation.
 *JCC(6/16/00)-add new status bits for tap ring
 *  U axis failure(ie. when the correction is applied) set bit 30 to 1(ON)
 *  V axis failure(ie. when the correction is applied) set bit 31 to 1(ON)
 *JCC(3/9/01)-update the document on 27th status bit which is now used
 *            by the tool 'hrc_correct_times' to indicate whether the time
 *            is corrected or not.
 *JCC(2/2002)-pass geompar to the pixlib call.
 *JCC(8/2002)-add amp_sf corrections; ( ampsfcorfile + do_amp_sf_cor ;);
 *JCC(11/2002)-add ra_nom/dec_nom; NOM_NOT_FOUND;
 *JCC(3/2003)-replace old hdrlib with new one 
 *            remove event_hdr.h/prim_hdr.h ; 
 *            remove PrimKey *primHdr_p; 
 *            replace EventKey* with Header_Type* for obs_info_p ;
 *JCC(4/2003)-new tapring spec ; 
 *           -read width_threshold from obs.par or evt.fits;
 *JCC(7/2003)-add a new function initial_status;
 *(4/2005)-comment out unused define's.
*(7/2005)-add removePath, stkExpand and ASOLFILE ;
*
*(5/2007)-new caldb4.
*(1/2009)-new gain 3dim image for hrcS ;  (obsolete 10/2009)
*        - add hpePrintErr();
*10/2009 - replace peterR hrcS 3dim gain image w/ dph hrcS gain table.
*10/2009 - add fap new hrcI gain image
*JCC(8/2012) - make TIMEGRID_LEN, RAWX_LEN dynamic for hrcS t_gain_map.
*    ( Note: the old 'fixed' values were TIMEGRID_LEN=18, RAWX_LEN=48 )
*************************************************************************/
#ifndef HRC_PROCESS_EVENTS_H
#define HRC_PROCESS_EVENTS_H

#include <stdio.h>
#include <bool.h>
#include <ctype.h>
#include <float.h>

#ifndef LIBERR_H
#include "liberr.h"
#endif

#ifndef STKLIB
#include "stack.h"
#define STKLIB
#endif

#ifndef PIXLIB_H
#include "pixlib.h"
#define PIXLIB_H
#endif

#ifndef ASCDM_H
#include "ascdm.h"
#define ASCDM_H
#endif

#ifndef HDRLIB2_H
#include "hdrlib2.h"
#endif

#include "badpixel_defs.h"

#ifndef ADCFILTER_DEFS_H
#include "adc_filter_defs.h"
#endif

#ifndef TAP_RING_DEFS_H
#include "tap_ring_defs.h"
#endif

#ifndef AMP_SF_COR_H
#include "amp_sf_cor_defs.h"
#endif

#ifndef L1_ALIGNMENT_DEFS_H
#include "l1_alignment_defs.h"
#define L1_ALIGNMENT_DEFS_H
#endif 

#ifndef L1_ASPECT_DEFS_H
#include "l1_aspect_defs.h"
#define L1_ASPECT_DEFS_H
#endif 

#ifndef L1H_DEGAP_DEFS_H
#include "l1h_degap_defs.h"
#define L1H_DEGAP_DEFS_H
#endif

#include "caldb4.h"

/* Here you must convert the name of your main function into an SPP/fortran
 * converted name.  Here is the rule.  Take the first 5 characters of 
 * your function name (that are literal characters, underscores do not count!)
 * and the last character, and follow it by an underscore.
 */
#define hrc_process_events evthrt_



/*  the following defines are used to set various debugging levels for the
 *  hrc_process_events module. These level values are compared against the 
 *  debug parameter of hrc_process_events.par file to determine the extent of 
 *  debugging information to log to a 'runlog' file.
 *
 *  DEBUG LEVELS  
 */
 
#define DEBUG_LEVEL_0     0   /* default - do not log any debugging info */ 
#define DEBUG_LEVEL_1     1   /* only log major info                     */  
#define DEBUG_LEVEL_2     2   /* log input and output attributes         */  
#define DEBUG_LEVEL_3     3   /* log statistical counts                  */ 
#define DEBUG_LEVEL_4     4   /* log almost everything                   */  
#define DEBUG_LEVEL_5     5   /* log everything                          */  


#define HDET_LOG_TO_STDOUT "stdout" /* log goes to screen instead of a file */
#define NEG_9999        -9999
#define NOM_NOT_FOUND   -9999.0    /* ra/dec_nom not found in obsfile */

/*----------- 10/2009 + 8/2012 - dph hrcS gain table for hrcS */
#define ORDER_LEN  2      /* polynomial order dim kk */
#define RAWY_LEN  576     /* RAWYGRID dim */
/* marcro: array index from n-dim to 1dim.
 *    gmIdx: for gainmap[kk=0:1, xj=0:47, yi=0:575];
 *    tgIdx: for tgain[yi=0:575, ww=0:17]
 *    g2Idx: for G_2nd[xj=0:47,  yi=0:575] 
 */
#define tgIdx(yi,ww)    ( (yi) + (ww)*RAWY_LEN )
#define gmIdx(kk,xj,yi,RAWX_LEN) ((kk)+(xj)*ORDER_LEN+(yi)*ORDER_LEN*RAWX_LEN)
#define g2Idx(xj,yi,RAWX_LEN)    ((xj)+(yi)*RAWX_LEN)
/*----------- end:*/


/*10/2009 - for inp_p->gainflag*/
#define OLD_SI_GAIN  0      /* old 2dim gain image for hrcS+I */
#define NEW_S_GAIN   1      /* new dph hrcS gain table */
#define NEW_I_GAIN   2      /* new fap hrcI gain image */

/*  the following is a list of all of the possible fields which may exist 
 *  in the level 1 event file. These defines are used to map input event
 *  file columns to an event data stucture which will hold all data 
 *  required for detector coordinate calculations on given events. The
 *  defines will also be used for mapping the data from the event data
 *  structure to the output event file.
 *
 *  HRC EVENT FILE COLUMNS
 */   

#define HDET_TIME               0  /* time of event                     */
#define HDET_MJR_FRAME          1  /* major frame                       */ 
#define HDET_MNR_FRAME          2  /* minor frame                       */
#define HDET_EVENT              3  /*                                   */  
#define HDET_CP_X               4  /* coarse x position                 */  
#define HDET_CP_Y               5  /* coarse y position                 */  
#define HDET_AX_1               6  /* x plane amp 1                     */  
#define HDET_AX_2               7  /* x plane amp 2                     */  
#define HDET_AX_3               8  /* x plane amp 3                     */  
#define HDET_AY_1               9  /* y plane amp 1                     */  
#define HDET_AY_2              10  /* y plane amp 2                     */  
#define HDET_AY_3              11  /* y plane amp 3                     */  
#define HDET_PHA               12  /*                                   */  
#define HDET_VETO_STATUS       13  /*                                   */  
#define HDET_EVENT_STATUS      14  /*                                   */  
#define HDET_X_POS             15  /*                                   */  
#define HDET_Y_POS             16  /*                                   */  
#define HDET_SUMAMPS           17  /*                                   */  
#define HDET_CHIP_ID           18  /*                                   */  
#define HDET_DUMMY             19  /*                                   */  
#define HDET_PI                20  /* pulse invarience                  */  
#define HDET_PHASCALE          21  /*                                   */  
#define HDET_RAWPHA            22  /*                                   */  
#define HDET_RAW_X             23  /* raw x coordinate                  */  
#define HDET_RAW_Y             24  /* raw y coordinate                  */  
#define HDET_DET_X             25  /* detector x coordinate             */  
#define HDET_DET_Y             26  /* detector y coordinate             */  
#define HDET_CHIP_X            27  /* chip x coordinate                 */  
#define HDET_CHIP_Y            28  /* chip y coordinate                 */  
/* #define HDET_AMP_TOT_X  29 */ /* sum of amps in x plane            */ 
/* #define HDET_AMP_TOT_Y  30 */ /* sum of amps in y plane            */ 
#define HDET_TICK              31  /*                                   */
#define HDET_SCIFR             32  /*                                   */
#define HDET_EVTCTR            33  /*                                   */
#define HDET_UNKNOWN_FIELD     34  /* unrecognized field                */ 
#define HDET_TDET_X            35  /* tdet x coordinate                 */
#define HDET_TDET_Y            36  /* tdet y coordinate                 */
#define HDET_SKY_X             37  /* sky (xy) x coordinate             */ 
#define HDET_SKY_Y             38  /* sky (xy) y coordinate             */ 
#define HDET_TAN_X             39  /* x axis tangent plane coordinate   */
#define HDET_TAN_Y             40  /* y axis tangent plane coordinate   */
#define HDET_FPZ               41  /* focal plane coord z axis          */  
#define HDET_STATUS            42  /* event status mask                 */ 
#define HDET_AMP_SF            43  /* scale used in degap/ratio checks  */
#define HDET_E_TRIG            44  /* combined u_trig and v_trig        */
#define HDET_VETO_STT          45  /* flight vstatus                    */
#define HDET_SUBMJF            46  /* sub major frame                   */
#define HDET_DET_ID            47  /* 0 = imaging, 1 = spectroscopy     */
#define HDET_MRF               48  /* */
#define HDET_STOPMNF           49 

#define HDET_SKY               50  /* vector column- x,y                */ 
#define HDET_DET               51  /* vector column- detx,dety         */ 
#define HDET_TDET              52  /* vector column- tdetx,tdety        */ 
#define HDET_CHIP              53  /* vector column- chipx,chipy        */ 
#define HDET_RAW               54  /* vector column- rawx,rawy          */ 

/* #define HDET_NUM_COLS  55 */ /* total number of possible fields   */



/*  the following defines represent the possible string column names that
 *  will be used to map event data to the appropriate fields in the event
 *  record structure 
 *
 *  VALID COLUMN NAMES
 */

#define HDET_TIME_COL          "TIME" 
#define HDET_MJR_FRAME_COL     "MJRFRAME" 
#define HDET_MJR_FRAME_COL_ALT1 "MJF" 
#define HDET_MNR_FRAME_COL     "MNRFRAME" 
#define HDET_MNR_FRAME_COL_ALT1 "MNF" 
#define HDET_MNR_FRAME_COL_ALT2 "STARTMNF" 
#define HDET_MNR_FRAME_COL_ALT3 "BEGINMNF" 
#define HDET_EVENT_COL         "EVENT"
#define HDET_CP_X_COL          "CP_X" 
#define HDET_CP_X_COL_ALT1     "CRSU"      /* crsu: coarse x coord. */
#define HDET_CP_X_COL_ALT2     "COARSEX"
#define HDET_CP_X_COL_ALT3     "COARSE_X" 
#define HDET_CP_Y_COL          "CP_Y" 
#define HDET_CP_Y_COL_ALT1     "CRSV"      /* crsv: coarse y coord. */
#define HDET_CP_Y_COL_ALT2     "COARSEY"
#define HDET_CP_Y_COL_ALT3     "COARSE_Y" 
#define HDET_AX_1_COL          "AU1" 
#define HDET_AX_1_COL_ALT1     "AX1" 
#define HDET_AX_2_COL          "AU2" 
#define HDET_AX_2_COL_ALT1     "AX2" 
#define HDET_AX_3_COL          "AU3" 
#define HDET_AX_3_COL_ALT1     "AX3" 
#define HDET_AY_1_COL          "AV1" 
#define HDET_AY_1_COL_ALT1     "AY1" 
#define HDET_AY_2_COL          "AV2" 
#define HDET_AY_2_COL_ALT1     "AY2" 
#define HDET_AY_3_COL          "AV3" 
#define HDET_AY_3_COL_ALT1     "AY3" 
#define HDET_PHA_COL           "PHA"
#define HDET_V_STS_COL         "VETOSTAT"
#define HDET_V_STS_COL_ALT1    "VSTAT"
#define HDET_E_STS_COL         "ESTAT"
#define HDET_STATUS_COL        "STATUS"
#define HDET_STATUS_COL_ALT1   "STAT"
#define HDET_X_POS_COL         "X_POS"
#define HDET_Y_POS_COL         "Y_POS"
#define HDET_SUMAMPS_COL       "SUMAMPS"
#define HDET_CHIP_ID_COL       "CHIP_ID"
#define HDET_DUMMY_COL         "DUMMY"
#define HDET_PI_COL            "PI"
#define HDET_PHASCALE_COL      "PHASCALE"
#define HDET_RAWPHA_COL        "RAW_PHA"
#define HDET_RAW_X_COL         "RAWX" 
#define HDET_RAW_X_COL_ALT1    "RAW_X" 
#define HDET_RAW_Y_COL         "RAWY" 
#define HDET_RAW_Y_COL_ALT1    "RAW_Y" 
#define HDET_DET_X_COL         "DETX" 
#define HDET_DET_X_COL_ALT1    "DET_X" 
#define HDET_DET_Y_COL         "DETY" 
#define HDET_DET_Y_COL_ALT1    "DET_Y" 
#define HDET_CHIP_X_COL        "CHIPX" 
#define HDET_CHIP_X_COL_ALT1   "CHIP_X" 
#define HDET_CHIP_Y_COL        "CHIPY" 
#define HDET_CHIP_Y_COL_ALT1   "CHIP_Y" 
#define HDET_TDET_X_COL        "TDETX"
#define HDET_TDET_X_COL_ALT1   "TDET_X" 
#define HDET_TDET_Y_COL        "TDETY"
#define HDET_TDET_Y_COL_ALT1   "TDET_Y" 
/* #define HDET_TAN_X_COL         "TANX"  */
/* #define HDET_TAN_X_COL_ALT1    "TAN_X" */
/* #define HDET_TAN_Y_COL         "TANY"  */
/* #define HDET_TAN_Y_COL_ALT1    "TAN_Y" */
#define HDET_SKY_X_COL         "X"
#define HDET_SKY_Y_COL         "Y"
#define HDET_TICK_COL          "TICK"
#define HDET_TICK_COL_ALT1     "CLKTICKS"
#define HDET_SCIFR_COL         "SCIFR"
#define HDET_EVTCTR_COL        "EVTCTR"
#define HDET_FPZ_COL           "FPZ"
#define HDET_AMP_SF_COL        "AMP_SF" 
#define HDET_E_TRIG_COL        "E_TRIG"
#define HDET_E_TRIG_COL_ALT1   "ETRIG"
#define HDET_VETO_STT_COL      "VETOSTT"
#define HDET_VETO_STT_COL_ALT1 "VETO_STT"
#define HDET_SUBMJF_COL        "SUBMJF"
#define HDET_SUBMJF_COL_ALT1   "SUB_MJF"
#define HDET_DET_ID_COL        "DET_ID"
#define HDET_MRF_COL           "MRF"
#define HDET_STOPMNF_COL       "STOPMNF" 
#define HDET_STOPMNF_COL_ALT1  "ENDMNF" 

#define HDET_SKY_COL           "SKY"
#define HDET_DET_COL           "DET"
#define HDET_TDET_COL          "TDET"
#define HDET_CHIP_COL          "CHIP"
#define HDET_RAW_COL           "RAW"

/*  the following defines are used by the coordinate transformation functions
 *  in hrc_process_events to represent the various coordinate systems.
 *
 *  COORDINATE SYSTEMS
 */
 
#define HDET_NO_COORD     "none"
#define HDET_COARSE_COORD "coarse"
#define HDET_CHIP_COORD   "chip"
#define HDET_TDET_COORD   "tdet"
#define HDET_DET_COORD    "det"
#define HDET_TAN_COORD    "tan"
#define HDET_SKY_COORD    "sky"
 
#define HDET_INV_VAL      -1
#define HDET_NONE_VAL      0
#define HDET_COARSE_VAL    1 
#define HDET_CHIP_VAL      2 
#define HDET_TDET_VAL      4 
#define HDET_DET_VAL       8 
#define HDET_TAN_VAL      16 
#define HDET_SKY_VAL      32 
 

/*  the following macros define the ranges of various arrays used in the
 *  hrc_process_events code. They are meant to assist in the readability of 
 *  the code. If these values are changed to adjust structure sizes, the 
 *  corresponding code must be compensated as well.
 *
 *  MISC. DEFINES 
 */ 

#define HDET_PLANE_X            0  /* x axis                             */ 
#define HDET_PLANE_Y            1  /* y axis                             */ 
#define HDET_PLANE_Z            2  /* z axis                             */
#define HDET_NUM_PLANES         2  /* total number of axes               */ 
#define HDET_1ST_AMP            0  /* first amp for a given axis         */ 
#define HDET_2ND_AMP            1  /* second amp for a given axis        */ 
#define HDET_3RD_AMP            2  /* third amp for a given axis         */ 
#define HDET_NUM_AMPS           3  /* number of amps per axis            */
#define HDET_1ST_ORD_CF         0  /* first order correction factor      */
#define HDET_2ND_ORD_CF         1  /* second order correction factor     */
#define HDET_NUM_CORR_FACTORS   2  /* number of correction factors       */ 
/* #define HDET_COL_DEF_LEN  256 */ /* length of evt column list string  */
#define TIME_INT_EPSILON        1  /* epsilon for tstart/tstop range     */ 
/* #define HDET_NUM_TAPS   256 */ /* number of tap positions  */
#define HDET_MAX_PI_NEW 1023 /*10/2009-max outCol PI for new hrcS + hrcI gain file*/
#define HDET_MAX_PI_OLD    255    /* max outCol PI for old gain img */

#define HRC_DEGAP_FILE_LEN  DS_SZ_PATHNAME /*length of the degap data file*/ 

#define GAIN_DEFAULT_CDELT_X  128  /* default for gain image (x axis)    */
#define GAIN_DEFAULT_CDELT_Y  128  /* default for gain image (y axis)    */

/* ----- new 3dim gain image */
#define   NEW_GAIN_TAPSIZE          256
#define   NEW_GAIN_SUBTAPS            3

#define HDET_WIRE_ON            0  /* wire charge test on                */ 
/* #define HDET_WIRE_OFF -1 */ /* wire charge test off  */ 

#define HRC_IMAGE_INST          0  /* hrc imaging instr. used */
/* #define HRC_SPECT_INST 1 */ /* hrc spectral instr. used */ 

/* #define HRC_SYS_LEN 10  */
/* #define BLANK_PARM " " */  /* 'empty' string to pass into pixlib */


/*  #define X_AXIS_KEY    "XS-INDXX"  */
/*  #define Y_AXIS_KEY    "XS-INDXY"  */
#define TELESCOP_KEY  "TELESCOP"
/* #define AC_SYSTEM     "ACSYS" */
#define DATAMODE_KEY  "DATAMODE"

/* #define NOM_RA_KEY    "RA_NOM" */
/* #define NOM_DEC_KEY   "DEC_NOM" */
/* #define NOM_ROLL_KEY  "ROLL_NOM" */

#define HRC_PROC_FLIGHT       0    /* process data as flight type */ 
#define HRC_PROC_XRCF         1    /* process data as xrcf type */ 

#define TELESCOP_XRCF_KEY   "XRCF-HRMA" 
#define TELESCOP_XRCF_KEY2   "XRCF/HRMA" 
#define NEXT_IN_LINE_VAL     "NEXT_IN_LINE" 
#define NEXT_IN_LINE_VAL2    "NEXT-IN-LINE" 

/*  the following defines are used by hrc_process_events to perform staging corrections
 *  on hsi data.
 *
 *  STAGING MACROS
 */

/* #define HSI_INSTRUMENT "hsi" *//* keyword to identify hsi instrument   */
#define STG_X_KEY          "STF_X"  /* keyword to identify x axis stage val */
#define STG_Y_KEY          "STF_Y"  /* keyword to identify y axis stage val */
#define STG_Z_KEY          "STF_Z"  /* keyword to identify z axis stage val */
#define STG_AIMPOINT_X           0  /* index of staging data x axis         */ 
#define STG_AIMPOINT_Y           1  /* index of staging data y axis         */ 
#define STG_AIMPOINT_Z           2  /* index of staging data z axis         */ 
#define HPY_YAW_KEY     "HRMA_YAW"  /* keyword for HRMA pitch yaw           */
#define HPY_PIT_KEY     "HRMA_PIT"  /* keyword for HRMA pitch yaw           */
#define HPY_PIT                  0  /* index of hpy pit data                */
#define HPY_YAW                  1  /* index of hpy yaw data                */

#define STG_ANG1_KEY    "STF_ANG1"  /* keyword to identify fam az           */ 
#define STG_ANG2_KEY    "STF_ANG2"  /* keyword to identify fam el           */ 
#define STG_ANG3_KEY    "STF_ANG3"  /* keyword to identify fam roll         */ 
#define STG_ANG1_NDX             0  /* index of fam az                      */ 
#define STG_ANG2_NDX             1  /* index of fam el                      */ 
#define STG_ANG3_NDX             2  /* index of fam roll                    */ 

#define SIM_X_KEY          "SIM_X"  /* keyword to identify sim x position   */ 
#define SIM_Y_KEY          "SIM_Y"  /* keyword to identify sim x position   */ 
#define SIM_Z_KEY          "SIM_Z"  /* keyword to identify sim x position   */ 
#define SIM_X_NDX                0  /* index of sim x pos                   */
#define SIM_Y_NDX                1  /* index of sim y pos                   */
#define SIM_Z_NDX                2  /* index of sim z pos                   */

/*  the following macros define the values of the various ARD file
 *  principal extension names. The values are from the HRC level 1 
 *  Pipeline ARD rev 1.4  02-22-99. 
 *
 *  ARD EXTENSION NAMES
 */

#define HPE_DEGAP_EXTNAME    "AXAF_DEGAP"
/* #define HPE_GAIN_EXTNAME     "AXAF_GAINMAP" */
#define HPE_ADC_EXTNAME      "AXAF_ADC"
#define HPE_HYP_EXTNAME	     "AXAF_FPTEST"	
#define HPE_AMP_FLAT_EXTNAME "AXAF_EFTEST"
#define HPE_AMP_SAT_EXTNAME  "AXAF_SATTEST"
 


/*  the following defines are used as bit masks to determine if all data
 *  dependencies are met in order to do detector coordinate calculations.
 *  the masks are all cast to unsigned short since that is the type of the
 *  actual variable being masked. An unsigned char probably would have
 *  sufficed for the current number of dependencies which exist, but the
 *  unsigned short was used to allow for future expansion.  
 *
 *
 *      bit  15 14 13 12 11 10  9  8  7  6  5  4  3  2  1  0
 *          |__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
 *                     \  \  \  \  \  \  \  \  \  \  \  \  \
 *                      \  \  \  \  \  \  \  \  \  \  \  \  \_ coarse x
 *           chipid     _\  \  \  \  \  \  \  \  \  \  \  \___ coarse y
 *           tdet y     _____\  \  \  \  \  \  \  \  \  \_____ amp 1 x axis
 *           tdet x     _________\  \  \  \  \  \  \  \_______ amp 2 x axis
 *           chip y     _____________\  \  \  \  \  \_________ amp 3 x axis
 *           chip x     _________________\  \  \  \___________ amp 1 y axis
 *                                           \  \_____________ amp 2 y axis
 *                                            \_______________ amp 3 y axis
 * 
 *  DATA DEPENDENCY MASK
 */

#define HDET_MASK_INIT   ((unsigned short) 0x00)   /* no dependencies met */
#define HDET_MASK_CPX    ((unsigned short) 0x01)   /* coarse x dep. met   */  
#define HDET_MASK_CPY    ((unsigned short) 0x02)   /* coarse y dep. met   */  
#define HDET_MASK_AX1    ((unsigned short) 0x04)   /* amp x 1  dep. met   */  
#define HDET_MASK_AX2    ((unsigned short) 0x08)   /* amp x 2  dep. met   */  
#define HDET_MASK_AX3    ((unsigned short) 0x10)   /* amp x 3  dep. met   */  
#define HDET_MASK_AY1    ((unsigned short) 0x20)   /* amp y 1  dep. met   */  
#define HDET_MASK_AY2    ((unsigned short) 0x40)   /* amp y 2  dep. met   */  
#define HDET_MASK_AY3    ((unsigned short) 0x80)   /* amp y 3  dep. met   */  
#define HDET_MASK_CHIPX  ((unsigned short) 0x100)  /* chip x column       */
#define HDET_MASK_CHIPY  ((unsigned short) 0x200)  /* chip y column       */
#define HDET_MASK_TDETX  ((unsigned short) 0x400)  /* tdet x column       */
#define HDET_MASK_TDETY  ((unsigned short) 0x800)  /* tdet y column       */
#define HDET_MASK_CHIPID ((unsigned short) 0x1000) /* chip id column      */

#define HDET_MASK_GOOD   ((unsigned short) 0xFF)    /* tdet dependency mask  */
#define HDET_MASK_CRS_REQ  ((unsigned short) 0xFF)  /* coarse dependency mask*/
#define HDET_MASK_CHIP_REQ ((unsigned short) 0x300) /* chip dependency mask  */
#define HDET_MASK_SCHIP_REQ ((unsigned short) 0x1300) /* chip dependency msk */
#define HDET_MASK_TDET_REQ ((unsigned short) 0xC00) /* tdet dependency mask  */
 


/*  The following macro's are used to keep track of the various output
 *  coordinate columns and whether to write them out as shorts or ints.
 *
 *  OUTPUT COORD COLUMN TYPES 
 */

/* #define HDET_COL_SKY_X               0  */
/* #define HDET_COL_SKY_Y               1  */
/* #define HDET_COL_DET_X               2  */
/* #define HDET_COL_DET_Y               3  */
/* #define HDET_COL_TDET_X              4  */
/* #define HDET_COL_TDET_Y              5  */
#define HDET_NUM_OUTCOORDS              6  

/* #define HDET_SHORT_COL               2  */
/* #define HDET_LONG_COL                4  */



/*  the following macro's define the bits used for the status mask of events
 *  processed by hrc_process_events. 
 *
 *     |3|3|2|2|2|2|2|2|2|2|2|2|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|0|0|0|0|0| 
 *     |1|0|9|8|7|6|5|4|3|2|1|0|9|8|7|6|5|4|3|2|1|0|9|8|7|6|5|4|3|2|1|0|
 * 
 *     +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+   
 *     | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |   
 *     +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+   
 *
 *     |   |     | | E |  VSTAT        | Process Stats |   | Data Stats| 
 *       \   \
 *        \   \
 *         \   \_ unused
 *          \
 *           \_ Tap ring status
 *
 *
 *
 *
 *
 *     31   Tap ring test failed on V
 *     30   Tap ring test failed on U
 *     29 - 28       - unused 
 *     27   Time correction status (1: time is changed ; 0: no change) 
 *     26   events telemetered in NEXT-IN-LINE mode 
 *     25   E status - V axis not triggered  
 *     24   E status - U axis not triggered
 *     23   V status - V axis center blank event 
 *     22   V status - U axis center blank event 
 *     21   V status - V axis width exceeded  
 *     20   V status - U axis width exceeded  
 *     19   V status - shield PMT active  
 *     18   V status - spare  
 *     17   V status - upper level discriminator exceeded  
 *     16   V status - lower level discriminator not exceeded  
 *     15   P status - Hot spots / events in bad region  
 *     14   P status - fine pos  
 *     13   P status - V center (bad amps)  
 *     12   P status - U center (bad amps) 
 *     11   P status - pha ratio  
 *     10   P status - zero psum  
 *     09   P status - grid ratio    
 *     08   P status - zero sum  
 *     07   PI > 255||1023     (old gain image|| new gain file)(10/2009)
 *     06   Out of sequence event 
 *     05   Amp flatness test failed on V
 *     04   Amp flatness test failed on U
 *     03   Amp saturation test failed on V
 *     02   Amp saturation test failed on U
 *     01   Hyperbolic test failed on V
 *     00   Hyperbolic test failed on U
 *
 *  STATUS MASK
 */

typedef long HRC_STATUS_T; 

#define HDET_V_TAP_RING_STS    ((HRC_STATUS_T) 0x80000000)  /*JCC: 2**31 */
#define HDET_U_TAP_RING_STS    ((HRC_STATUS_T) 0x40000000)  /*JCC: 2**30 */

#define HDET_NEXT_IN_LINE_STS  ((HRC_STATUS_T) 0x04000000) 

/* #define HDET_V_TRGGR_STS       ((HRC_STATUS_T) 0x02000000)  */
/* #define HDET_U_TRGGR_STS       ((HRC_STATUS_T) 0x01000000)  */

/* #define HDET_V_BLNK_CTR_STS    ((HRC_STATUS_T) 0x00800000)   */
/* #define HDET_U_BLNK_CTR_STS    ((HRC_STATUS_T) 0x00400000)   */

#define HDET_V_WIALIGN_STS     ((HRC_STATUS_T) 0x00200000) 
#define HDET_U_WIALIGN_STS     ((HRC_STATUS_T) 0x00100000) 
/* #define HDET_SHLD_PMT_STS      ((HRC_STATUS_T) 0x00080000) */
/* #define HDET_SPARE_STS         ((HRC_STATUS_T) 0x00040000)  */
/* #define HDET_UP_DSCRM_STS      ((HRC_STATUS_T) 0x00020000)  */
/* #define HDET_LO_DSCRM_STS      ((HRC_STATUS_T) 0x00010000)  */

#define HDET_HOT_SPOT_STS      ((HRC_STATUS_T) 0x00008000) 
#define HDET_FIN_POS_STS       ((HRC_STATUS_T) 0x00004000) 
#define HDET_V_CNTR_STS        ((HRC_STATUS_T) 0x00002000) 
#define HDET_U_CNTR_STS        ((HRC_STATUS_T) 0x00001000) 
#define HDET_PHA_RTO_STS       ((HRC_STATUS_T) 0x00000800) 
#define HDET_ZERO_PSUM_STS     ((HRC_STATUS_T) 0x00000400) 
#define HDET_GRID_RTO_STS      ((HRC_STATUS_T) 0x00000200) 
#define HDET_ZERO_SUM_STS      ((HRC_STATUS_T) 0x00000100) 
#define HDET_PI_VALUE_STS      ((HRC_STATUS_T) 0x00000080) 
#define HDET_SEQUENCE_STS      ((HRC_STATUS_T) 0x00000040) 

#define HDET_V_AMP_FLAT_STS    ((HRC_STATUS_T) 0x00000020) 
#define HDET_U_AMP_FLAT_STS    ((HRC_STATUS_T) 0x00000010) 
#define HDET_V_AMP_SAT_STS     ((HRC_STATUS_T) 0x00000008) 
#define HDET_U_AMP_SAT_STS     ((HRC_STATUS_T) 0x00000004) 
#define HDET_V_HYP_STS         ((HRC_STATUS_T) 0x00000002) 
#define HDET_U_HYP_STS         ((HRC_STATUS_T) 0x00000001) 



/*  the following structure holds all information pertaining to any single
 *  event being processed. Data from the input file is loaded into this
 *  structure and is processed. Relevant calculations such as raw x and
 *  y coordinates are stored after being calculated. Upon completion of 
 *  processing, hrc_process_events writes out specified fields from this structure
 *  to an output event file. These output fields are cast as appropriate.
 *
 *  EVENT RECORD STRUCTURE
 */

typedef struct event_rec_t {
   double time;                               /* time of event           */ 
   double amps_dd[HDET_NUM_PLANES][HDET_NUM_AMPS];/*amplitudes for computation*/
   HRC_STATUS_T status;                       /* event status mask       */ 
   long    tick;                              /* sac                     */ 
   long    scifr;                             /* sac                     */ 
   long    major_frame;                       /* major frame             */
   long    minor_frame;                       /* minor frame             */
   long gain_index[2]; /*used for old 2dim gain image; not for dph hrcS gain table*/
   short   stopmnf; 
   short   mrf;   
   short  submjf;                             /* sub major frame         */
   short  event;
   short  cp[HDET_NUM_PLANES];     /* crsu,crsv: x,y coarse positions*/ 

  /**********************************************************************
   * JCC(5/11/00)
   *
   * short amps_sh- original raw amplitudes from infile and will be in outfile.
   *                ie.  (AU1, AU2, AU3) are raw amplitudes for X/U plane 
   *                     (AV1, AV2, AV3) are raw amplitudes for Y/V plane 
   *
   * double amps_dd - modified amplitudes (eg. AU3_tap_corr); used for computation; 
   *
   * short sum_amps = io_col-sumamps = AU1+ AU2+ AU3_tapCorr+ AV1+ AV2+ AV3_tapCorr
   **********************************************************************/
   short  amps_sh[HDET_NUM_PLANES][HDET_NUM_AMPS]; 
   /*short amps_3RD_raw[HDET_NUM_PLANES]; *//*original 3RD amps: AU3 & AV3 */
   short  amps_tap_flag[HDET_NUM_PLANES]; /* applying 'tap correction' ? */
   short  pha;
   unsigned short  event_status;
   short  xpos;
   short  ypos;
   unsigned short  sum_amps;                  /* io_column 'sumamps' */ 
   short  chipid;                             /* identifies which chip   */
   short  dummy;
   long   pi;        /*(1/2009-outCol PI; pulse invarience*/ 
double pi_double; /*10/2009-interm. var. to calc. pi for dph hrcS gainTab; spi;*/
double DDn;       /*10/2009 - to compute pi_double for dph hrcS gain table*/
   short  phascale;
   short  rawpha;    /* RAWX, RAWY for the output */
   short  evtctr;                             /* sac                     */ 
   short  amp_sf;      /* output is the modified amp_sf; "scale" used in degap*/ 
   VEC2_DBLE amp_tot;      /* sum of amps per plane after tap correctins */ 
   VEC2_DBLE rawpos;       /* io raw x and y coordinates; outCol may change;*/ 
   VEC2_DBLE tdetpos;                         /* tdetector x and y coords*/ 
   VEC2_DBLE detpos;                          /* detector x and y coords */ 
   VEC2_DBLE chippos;                         /* chip x and y coordinates*/ 
   VEC2_DBLE workpos;                         /* radomized chip coords   */ 
   VEC2_DBLE skypos;                          /* "sky" coords            */ 
   VEC2_DBLE fppos;                           /* focal plane coordinates */ 
   VECS_CHAR chipname;   
   unsigned char veto_stt;                    /* veto status (flight)    */
   unsigned char veto_status;                 /* veto status (lab/hsi)   */
   unsigned char e_trig;                      /* trigger (flight)        */ 
   unsigned char det_id;                      /* 0 = imaging 1 = spect   */ 

} EVENT_REC_T, *EVENT_REC_P_T;


/*  The following structure is used by hrc_process_events to store information
 *  pertaining to event files (both input and output). The fields in this
 *  structure are populated via a call to setup_event_files().
 *
 *  EVENT FILE SETUP STRUCTURE
 */
typedef struct event_setup_t {
   dmDataset*    dataset;    /* dataset (file) handle                      */
   dmBlock*      primary;    /* header keyword extension handle            */
   dmBlock*      extension;  /* event data extension handle                */
   dmDescriptor**  desc;     /* event attributes                           */
   long*         dim;        /* column dimension (> 1 for vectors)         */
   int           num_cols;   /* number of columns in event file            */
   int           num_rows;   /* number of rows in the event file           */
   dmDataType*   types;      /* ptr to array of event data column types    */
   short*        mapping;    /* mapping of data structure and file columns */
   char*         file;       /* file name                                  */ 
   char*         eventdef;   /* output columns                             */
} EVENT_SETUP_T, *EVENT_SETUP_P_T;

  /*-----------------------------------------------------------
   * variable flg  :
   *    INIT_OK = callInit was called and is ok.
   *    INIT_NOT_OK = callInit was called but there's problem
   *    INIT_NOT_NEED = all input files are not 'caldb'.
   *----------------------------------------------------------*/
  typedef enum
  {
      INIT_OK=1,
      INIT_NOT_OK=2,     /* calInit error */
      INIT_NOT_NEED=3    /* no caldb interface */
  }   hrc_FLAG ;

  typedef struct hrc_for_caldb4
  {
          char   *telescop;
          char   *instrume;
          char   *tmp_detnam ;      /* hrc_si -> hrc_s */

      hrc_FLAG   flg ;              /* status of calling calInit */

      calCALDB   *myCaldb  ;
     calSEARCH   *mySearch ;
    Header_Type  *hdr ;

           int   debug ;

  } HRC_CALDB4_T,   *HRC_CALDB4_P ;

 
/*  the following structure is used to hold the input parameters which are
 *  used by hrc_process_events for various processing information. The values 
 *  contained in the structure are loaded in from the hrc_process_events 
 *  parameter file.
 *
 *  INPUT PARAMETERS STRUCTURE 
 */

typedef struct input_parms_t {
   double cf[HDET_NUM_PLANES][HDET_NUM_CORR_FACTORS];
   double obs_tstart;      /* tstart from obs.par file                       */
   double obs_tstop;       /* tstop from obs.par file                        */
   double evt_tstart;      /* start time (from earliest event)               */
   double evt_tstop;       /* end time (from latest   event)                 */
   double evt_mjd_obs;     /* mjd_obs key in evtfile */
   double pha_ratio;
   double grid_ratio;
   double amp_gain;        /* amp gain used for degap ratio checking         */
   double default_time;    /* default time field (primarily for HSI)         */
   double end_time;        /* end time for header (primarily for HSI)        */
   double time_offset;     /* time offset to synch to fam data times         */
   double fp_scale;        /* fp scale in arc seconds                        */
   double crval[3];        /* sky coordinate values for wcs/asp corrections  */
   double crpix[2];        /* sky coordinate values for wcs/asp corrections  */
   double cdelt[2];        /* sky coordinate values for wcs/asp corrections  */
   float  randpixsize;     /* width of randomization (-size..+size)          */
   unsigned long rand_seed; /* seed for pixlib randomization (0 = use time)  */
   long   gain_axlen[2];  /* axes lengths for old 2dim gain image */

/* ---- 10/2009 - for dph hrcS gain table. */
short gainflag; /*2=fap new hrcI gainImg; 1=dph hrcS gainTab; 0=old 2dim gainImg*/
long gainmapSize ;   /* gainmap column in new hrcS gainTab */
long tgainSize ;     /* tgain column in new hrcS gainTab */
long rawxgridSize ;  /* rawxgrid column in new hrcS gainTab*/
long rawygridSize ;  /* rawygrid column in new hrcS gainTab */
long timegridSize ;  /* timegrid column in new hrcS gainTab */
double *gainmapVal ;   /* gainmap column ; dim=2*48*576*/
double *tgainVal ;     /* tgain column ; dim=576*18*/ 
double *rawxgridVal ;  /* rawxgrid column ; dim=48 */
double *rawygridVal ;  /* rawygrid column ; dim=576 */
double *timegridVal ;  /* timegrid column ; dim=18 */
double *obs_tgain;     /* adjusted tgain(func. of evt MJD_OBS) ; dim=576*/ 
double *G_2nd;         /* adjusted gainmap ; dim=48*576 */
/* ---- end */

/*10/2009-value of key sampnorm in fap new hrcI gainImg (only); gain_factor;*/
double sampnorm;  

   int    debug;           /* specifies level of debug information to print  */
   int    wire_charge;
   short  coord_type[HDET_NUM_OUTCOORDS]; /* output col data type to utilize */
   boolean pix_init;       /* TRUE = pixlib has been configured              */
   boolean next_in_line;   /* TRUE = event data telemetered as next-in-line  */
   boolean do_raw;         /* TRUE = calculate raw coordinates               */
   boolean do_pi;          /* TRUE = calculate pulse invarience              */ 
   boolean clobber;        /* TRUE = remove existing output file w/ same name*/ 
   boolean scl_xsts;       /* TRUE = scale column (amp_sf) is in input file  */
   boolean do_ratio;       /* TRUE = perform ratio validity checks           */ 
   boolean do_ADC;         /* TRUE = perform ADC corrections                 */ 

   /* for amp_sf_cor */
   boolean do_amp_sf_cor;  /* TRUE = perform amp_sf correction (from hpe.par)*/
   boolean get_range_switch_level; /* TRUE=range_switch_level already got from
                                      obspar or event lev1 */
   boolean match_range_switch_level; /* TRUE=find the matched range_switch_level
                                        in amp_sf_corr caldb file */
   boolean AMPSFCOR;       /* TRUE=amp_sf_cor is applied to data ; */ 
   boolean evt_AMPSFCOR;   /* AMPSFCOR from evt1 infile; default=FALSE; */
   short  range_switch_level; /* range_switch_level from obs.par or 
                                 RANGELEV from evt1.fits*/
   short  obs_r_s_l ;   /* range_switch_level from obs.par; initial=-9999 */
   short  evt_r_s_l ;   /* RANGELEV from evt infile; initial=-9999 */

   short  obs_widthres ; /* -9999; width_threshold in obs.par==WIDTHRES evt1 */
   short  evt_widthres ;  /* width_threshold in obs.par = WIDTHRES evt fits */

   short  x_taps;          /* number of taps on x axis                       */ 
   short  y_taps;          /* number of taps on y axis                       */ 
   short  hrc_system;      /*0=image, 1=spect, 2=spect in imaging mode       */
   short  grating;         /* grating type - 0=none, 1=letg, 2=hetg, 3=toga  */ 
   short  start;           /* start of coordinate transformations            */
   short  stop;            /* end of coordinate transformations              */
   short  x_col;           /* column being used as x for coord calcs         */ 
   short  y_col;           /* column being used as x for coord calcs         */ 
   short  fpsys;           /* focal plane system ID                          */
   short  min_tap[HDET_NUM_PLANES]; /* min value for tap positions           */
   short  max_tap[HDET_NUM_PLANES]; /* max value for tap positions           */
   short  gain_cdelt[HDET_NUM_PLANES]; /* size of the gain image axes        */
   short  processing;      /* 1 = XRCF,   HRC_PROC_FLIGHT (0) = flight       */
   char*  old_sys;         /* hrc system used in previous run                */
   char*  ASOLFILE;        /* asol filenames without path */
   char   degap_file[DS_SZ_PATHNAME]; /* data to use for degap table     */ 
   char   stack_in[DS_SZ_PATHNAME]; /* I - file containing event file list   */
   char   geompar[DS_SZ_PATHNAME];  /* I - geompar for pixlib */
   char   outfile[DS_SZ_PATHNAME];  /* I - file name of output event file    */
   char   obsfile[DS_SZ_PATHNAME];  /* I - name of obs.par file              */
   char   logfile[DS_SZ_PATHNAME];  /* I - file name of output debug log file*/
   char   align_file[DS_SZ_PATHNAME]; /* I - path/name of alignment file     */
   char   asp_file[DS_SZ_PATHNAME]; /* I - path/name of aspect file          */
   char   gain_file[DS_SZ_PATHNAME]; /* I - path/name of gain image file     */
   char   adc_file[DS_SZ_PATHNAME]; /* I - path/name of adc correction file  */
   char   badpixfile[DS_SZ_PATHNAME]; /* I - path/name of bad pixel file     */
   char   hypfile[DS_SZ_PATHNAME]; /* I - path/name of hyperbolic test file  */
   char   tapfile[DS_SZ_PATHNAME]; /* I - path/name of tap ring test file*/
   char   ampsfcorfile[DS_SZ_PATHNAME]; /* I -path/name of amp_sf_correction */
   char   ampflatfile[DS_SZ_PATHNAME]; /* I - path/name of flatness test file*/
   char   ampsatfile[DS_SZ_PATHNAME]; /* I - path/name of saturation tst file*/
   char   time_start[DS_SZ_KEYWORD];  /* I - default start time value        */
   char   time_stop[DS_SZ_KEYWORD];  /* I - end time value for header        */
   char   outcols[DS_SZ_COMMAND];    /* I - list of columns to write out     */
   char   badfile[DS_SZ_PATHNAME];   /* I - name of output 'bad' event file  */
   char   badoutcols[DS_SZ_COMMAND];  /* I - columns to write to 'bad' file  */
   char   start_coord[DS_SZ_KEYWORD]; /* I - start of coord transformations  */
   char   stop_coord[DS_SZ_KEYWORD];  /* I - end point of coord transforms   */

   char   instrume[DS_SZ_KEYWORD];  /* instrument name- used as par file name*/
   char   datamode[DS_SZ_KEYWORD];  /* datamode keyword value              */
   char   telescop[DS_SZ_KEYWORD];   /* telescop keyword value             */
   Header_Type *obs_info_p;          /* pointer to the new header lib        */

   double   sim_x;                   /* sim_x value                          */
   double   sim_y;                   /* sim_y value                          */
   double   sim_z;                   /* sim_z value                          */
   double   ra_nom, dec_nom ;        /* crval[0..1] */
   short    mjd_obs_warn;     /*for dph hrcS gain table */
   short    use_obs;                 /* flag- 1 = obs read, 0 = no obs       */

   Header_Type   *caldb4_hdr;
   HRC_CALDB4_P  hcp ;

} INPUT_PARMS_T, *INPUT_PARMS_P_T;




/*  the following structure is used by hrc_process_events to keep track of 
 *  processing statistics. This values are output so that the user can 
 *  determine the extent of accuracy in processing. 
 *
 *  STATISTICS STRUCTURE 
 */

typedef struct statistics_t {
   long total_events_in;
   long total_events_out;
   long num_files_in;
   long num_bad_files;
   long bad_grid_ratio;
   long bad_pha_ratio;
   long bad_dist[HDET_NUM_PLANES];
   long bad_bot;
   long fixed_mfinpos;
   long fixed_pfinpos;
   long sequence_err; 
   unsigned short dependencies; /* dependency mask */ 
} STATISTICS_T, *STATISTICS_P_T;



#define HPE_LEN_80     80 

/*macro: if X>Y, return 1; */
#define hpe_gt(X,Y)   (((X)- (Y)) > (fabs(X)*2.0*DBL_EPSILON) ? 1:0 )

/*macro: if X>=Y, return 1; */
#define hpe_ge(X,Y)   (((X)- (Y)) >= (fabs(X)*2.0*DBL_EPSILON) ? 1:0 )


/*  the following typedef is used to store the various parameter file keywords
 *  which are read in from the instrument parameter file and written out to
 *  the output event file.
 *
 *  INSTRUMENT KEYWORDS STRUCTURE
 */

typedef struct inst_keywords_t {
   VECL_CHAR fp_system;
   VECL_CHAR chip_system;
   VECL_CHAR tdet_system;
} INST_KEYWORDS_T, *INST_KEYWORDS_P_T; 




/*  the following is a list of function prototypes of the routines
 *  that were created explicitly for hrc_process_events.
 *
 *  FUNCTION PROTOTYPES 
 */  

/* routine to process hrc events and return coordinate information */ 
extern dsErrCode hrc_process_events(void);

/* routine to compute pulse invarience for a given hrc event */ 
extern void   calculate_pi_hrc(float*, INPUT_PARMS_P_T, EVENT_REC_P_T); 

/* routine to get pulse height totals for a given hrc event */ 
extern void   sum_phas_hrc(EVENT_REC_P_T); 

/* routine to get pulse height totals for a given hrc event */ 
extern void   ratio_checks_hrc(EVENT_REC_P_T, 
                               INPUT_PARMS_P_T, 
                               STATISTICS_P_T); 

/* routine to write a given hrc event to an output file */ 
extern void   write_hrc_events(EVENT_SETUP_P_T, 
                               EVENT_REC_P_T,
                               dsErrList*);  

/* routine to set up mapping of event columns */
extern bool   parse_hrc_evt_columns (char**, 
                                     int, 
                                     short*); 

/* routine to map input event columns to event record structure */
extern void   load_event_data(EVENT_SETUP_P_T, EVENT_REC_P_T);
extern void   initial_status(INPUT_PARMS_P_T,  EVENT_REC_P_T); 

/* routine to setup bit mask for data dependency check */
extern unsigned short   dependency_check_init (short*, 
                                               int);

/* routine to determine if data dependencies are all met */
extern void   dependency_check_hrc (INPUT_PARMS_P_T, 
                                    STATISTICS_P_T,
                                    dsErrList*);

/* routine to ensure that output coordinate cols are calculated */
extern boolean   output_coord_validate(short*, 
                                       int, 
                                       INPUT_PARMS_P_T);

/* routine to read input parameters */ 
extern void   load_input_parameters(INPUT_PARMS_P_T, 
                                    dsErrList*);

/* routine to write out instrument parameter file keywords */ 
extern void   write_instrume_params(dmBlock*, 
                                    INST_KEYWORDS_P_T);

/* routine to load instrument parameter keywords into data structure */
extern boolean   read_instrume_params(char*, 
                                      INST_KEYWORDS_P_T);

/* routine to set staging/fam values for pixlib */
extern void   set_up_mirror(dmBlock*, 
                            INPUT_PARMS_P_T, 
                            char*, 
                            ALIGNMENT_REC_P_T,
                            dsErrList*);

/* routine to compute chip coords amd chip ids */
extern void   hrc_tel_to_chip(double*, 
                              short, 
                              short*, 
                              double*);

/* routine to add a hot spot to the ordered bad pixel list */ 
extern void   update_bad_pixel_list(BAD_PIX_P_T, 
                                    BAD_PIX_P_T*);
 
/* routine to remove the hot spot list and deallocate the list's memory */ 
extern void   cleanup_bad_pixel_data(BAD_PIX_A_T);
 
/* routine to load a stack of hot spot (bad pixel) files into memory */ 
extern boolean   load_bad_pixel_files(char*, 
                                      BAD_PIX_P_T*);
 
/* routine to access and load a single hot spot (bad pixel) file */ 
extern boolean   open_bad_pixel_file(BAD_PIX_SETUP_P_T);
 
/* routine to parse a hot spot (bad pixel) column name */ 
extern short   map_bad_pixel_column(char[]);
 
/* routine to load a column value from a row in a hot spot file */
extern void   load_bad_pixel_data(BAD_PIX_SETUP_P_T, 
                                  BAD_PIX_P_T*);

/* function to deallocate memory and close descriptors for bad pixel routines*/
extern void   close_bad_pixel_file(BAD_PIX_SETUP_P_T hk_p);
 
/* function to check for hot pixels */
extern void   check_for_bad_pixels(BAD_PIX_A_T, 
                                   EVENT_REC_P_T);
 
/* function to adjust the output eventdef based on output system */
extern void   adjust_output_eventdef(INPUT_PARMS_P_T); 

extern void hpe_setup_degap_file(INPUT_PARMS_P_T,
                                 DEGAP_CONFIG_P_T*,
                                 dsErrList*);

/* function to associate internal defines to columns to use in transforms */
extern boolean   map_start_column(short, 
                                  short*, 
                                  short*); 

/* function to determine what coord systems to transform */ 
extern short   parse_coord_range(char*); 

/* routine to open the input file and take care of input housekeeping */
extern void   hrc_process_setup_input_file(EVENT_SETUP_P_T,
                                           INPUT_PARMS_P_T, 
                                           STATISTICS_P_T,
                                           dsErrList*);

/* routine to open the output file and take care of input housekeeping */
extern void   hrc_process_setup_output_file(EVENT_SETUP_P_T,
                                            EVENT_SETUP_P_T,
                                            INPUT_PARMS_P_T, 
                                            ALIGNMENT_REC_P_T, char***,
                                            dsErrList*);
extern void removePath(char *inStr, char **outStr);
extern void stkExpand( char *inStr, char **outStr );

/* function to free memory and close file descriptors */ 
extern void   hrc_process_evt_file_cleanup(EVENT_SETUP_P_T); 

/* routine to create output event columns and set column info keywords */
extern void hrc_setup_columns(EVENT_SETUP_P_T, INPUT_PARMS_P_T, char**, dsErrList*);

void hpeSetRange_s(dmBlock* bb, dmDescriptor* dd, char *n1, short v1, short v2);
void hpe_set_ranges( EVENT_SETUP_P_T evtout, INPUT_PARMS_P_T inp, char** name);

/* routine to load old or new gain file */
extern void load_gain_image(char*, INPUT_PARMS_P_T, float**, dsErrList*);
 
/* routine to setup debug/log file */
extern FILE* hrc_process_setup_logfile(INPUT_PARMS_P_T,
                                       dsErrList*);
 
/* routine to compute coordinate values for a given hrc event */
extern boolean calculate_coords_hrc(EVENT_REC_P_T,
                                    INPUT_PARMS_P_T,
                                    STATISTICS_P_T,
                                    ASPECT_REC_P_T,
                                    DEGAP_CONFIG_P_T,
				    short,
                                    dsErrList*);
 
/* routine to print out and remove warnings from error list */
extern boolean process_warnings(dsErrList*, FILE*, int);
 
/* routine to verify event times against obs.par tstart/tstop */
extern void hrc_process_time_check(INPUT_PARMS_P_T,
                                   FILE*,
                                   dsErrList*);

/* routine to read obs.par */
extern void hrc_process_read_obsfile(INPUT_PARMS_P_T,
                                     dsErrList*); 

/* routine to set the pixlib configuration for the current run */
extern void hrc_process_configure_pixlib(EVENT_SETUP_P_T,
                                         INPUT_PARMS_P_T,
                                         INST_KEYWORDS_P_T,
                                         ALIGNMENT_REC_P_T,
                                         dsErrList*); 

/* routine to set instrument specific values */
extern void hrc_process_set_instrume(INPUT_PARMS_P_T,
                                     dsErrList*);

extern void calc_coarse_coords(EVENT_REC_P_T,
                               INPUT_PARMS_P_T,
                               STATISTICS_P_T,
                               DEGAP_CONFIG_P_T,
                               dsErrList*);

boolean calc_chip_coords(EVENT_REC_P_T,
                         INPUT_PARMS_P_T,
                         STATISTICS_P_T,
                         DEGAP_CONFIG_P_T,
                         dsErrList*);

/* saturation test function */
extern void open_amp_saturation_file (INPUT_PARMS_P_T, 
                                      SAT_TEST_P_T *, 
                                      dsErrList *);

/* routine to open/read flatness test ARD */
extern void open_evt_flatness_file(char *, 
                                   double **, 
                                   dsErrList *);

/* routine to open/read hyperbolic test ARD */
extern void open_hyperbolic_file(char *, 
                                 HYP_TEST_P_T *, 
                                 dsErrList *);

/* routine to evaluate saturation test */
extern void check_amp_saturation(EVENT_REC_P_T, 
                                 SAT_TEST_P_T);

/* routine to evaluate flatness test */
extern void check_evt_flatness(EVENT_REC_P_T, 
                               double);
 
/* routine to evaluate hyperbolic test */
extern void check_hyperbolic(EVENT_REC_P_T, 
                             HYP_TEST_P_T);

/* helper function for the hyperbolic test */
extern double hyperbolic_help(double, 
                              double, 
                              double, 
                              double);

/* amp_sf correction functions */
extern void open_amp_sf_cor_file(INPUT_PARMS_P_T inp_p, 
             AMPSFCOR_COEFF_P_T *ampsfcor_coeff, dsErrList *hpe_err_p) ;

extern void sum_raw_amps(EVENT_REC_P_T evt_p, double *tot_raw_amps ) ; 

extern void apply_amp_sf_cor(EVENT_REC_P_T evt_p, 
             AMPSFCOR_COEFF_P_T ampsfcor_coeff ) ; 

extern void write_amp_sf_corr_key( dmBlock* extension, INPUT_PARMS_P_T inp_p);

/* tap correction functions */
extern void open_tap_ring_file(INPUT_PARMS_P_T, 
                  TRING_COEFFS_P_T *, dsErrList *);

extern void check_tap_ring(INPUT_PARMS_P_T,EVENT_REC_P_T, 
                           TRING_COEFFS_P_T ) ; 

/* for old 2dim gain image or fap new hrcI 2dim gain image*/
/* 10/2009 - rename old_gain_index to image_2dim_gain_index(); no code changes;*/
extern void image_2dim_gain_index( INPUT_PARMS_P_T inp_p, EVENT_REC_T *evt_p);

/* 10/2009 - functions for dph hrcS gain table */
extern void find_timegrid_val(INPUT_PARMS_P_T inp_p, long *ww, 
                              double *timeGrid,short *found);
extern long grid_idx( long dim_gainRaw, double *gainRaw, double evtRaw, 
                      short *found);
extern void calc_DDn ( EVENT_REC_T *evt_p) ;
extern void set_gainflag( INPUT_PARMS_P_T inp_p, dsErrList* err_p);
extern void read_gainTab_col(dmBlock* srcBlock, char* colName,
                             double** colVal, long* colSize);
extern void load_S_new_gain_table( INPUT_PARMS_P_T inp_p);
extern void S_new_gain_index_pi(INPUT_PARMS_P_T inp_p, EVENT_REC_T *evt_p);
/* end: */ 

/* 10/2009 - functions for both dph and fap new gain file*/
extern void check_spi_limit( EVENT_REC_P_T evt_p ) ;


/* display error msgs */
extern dsErrCode hpePrintErr( dsErrList* hpe_err_p, FILE* log_ptr, int debug);


/* Access the CALDB for all calibration files */
extern void access_caldb( INPUT_PARMS_P_T , dsErrList*);

   extern HRC_CALDB4_P  init_caldb_var ( INPUT_PARMS_P_T inp_p ) ;
   extern dsErrCode find_caldb_file( char *myFile, char *myProduct, HRC_CALDB4_P hcp);

#endif   /* closes #ifndef HRC_PROCESS_EVENTS_H */  
