/*                                                                
**  Copyright (C) 1997-2009,2010  Smithsonian Astrophysical Observatory 
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

/*H***********************************************************************
* FILE NAME: tap_ring_sub.c
* DEVELOPEMENT: tools
* DESCRIPTION: This file contains functions used for 'tap ring corrections'
*
* The purpose of tap corrections is to make correction for A3 and
* recompute the fine coordinates using the new A3. In addition, 
* the new A3 should be used in the other tests. But for the output 
* file, make sure to write the original A3.
*
* There will be no tap corrections if tapfile='NONE', or if it does
* not pass certain conditons described in 'check_tap_ring'.
*
* Notes :
*
*  'check_tap_ring()' returns a NEW amps_dd[][3RD_AMP]
*      evt_p->amps_dd[plane][1ST_AMP] 
*      evt_p->amps_dd[plane][2ND_AMP] 
*      evt_p->amps_dd[plane][3RD_AMP] = tap_corrected_A3
*  'sum_phas_hrc()' returns
*      evt_p->amp_tot[plane] = amps_dd[1ST] + amps_dd[2ND] + tap_corrected_A3'
*
* JCC(5/1/00) - initial version
* JCC(5/11/00)- use amps_dd for computation.
* JCC(6/16/00)-add new status bits for tap ring ( bad_bit_mask )
*  U axis failure(ie. when the correction is applied) set bit 30 to 1(ON)
*  V axis failure(ie. when the correction is applied) set bit 31 to 1(ON)
* JCC(2/9/01) -enhanced to memset axis_val_a  before calling dmGetScalar_c.
* JCC(6/7/01) -passing maxlen-1 to dmGetScalar_c
*JCC(4/2003)-new tapRing spec ;
*  - new column TRING_OFFSET in CALDB for hrcS and hrcI ;
*  - keyword width_threshold in obs.par == WIDTHRES in evt1 ; 
*    ( if(obs_widthres==evt_widthres=NEG_9999) -> ie. value missing in obs&evt
*        it will still go thru the condition check for tapRing corr. )
*JCC(6/2003)- if TRING_OFFSET column missing in CALDB, set value to 0.0 and
*     continue tapRing_corr.; result should be same as old tapRing spec.;
*JCC(7/2003)-reset status bit before reapply tap ring corr.;
* (11/2008) - add hpe_gt for bugfix.
(1/2009)-updated comments.
*H***********************************************************************/
#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif

/**********************************************************************
 * open_tap_ring_file :
 *    Get the column values from the tap ring coefficients file ; 
 *    There should be only 2 rows (U & V) in the file.  
 **********************************************************************/
void open_tap_ring_file ( 
        INPUT_PARMS_P_T   inp_p,            /*i: passing tapfile and debug*/
        TRING_COEFFS_P_T *tring_coeffs_p, /*o: 1 dim arr with only 1 elem*/
        dsErrList         *hpe_err_p        /*u */  )
{
   /* 1: get TRING_OFFSET value from column ; 0: set it to 0.0 */ 
   short set_TRING_OFFSET = 1 ; 

   *tring_coeffs_p = NULL;      /* a scalr holds 2 rows of data */

   /* make sure a file name is provided before processing */
   if(((strcmp(inp_p->tapfile, "NONE") != 0) &&
      (strcmp(inp_p->tapfile, "none") != 0)) &&
      (strcmp(inp_p->tapfile, "\0") != 0))
   {
      if (dmDatasetAccess(inp_p->tapfile, "R") == dmTRUE)
      {
         short ii;
         short max_rows = 2;
         char  col_names_a[TRING_COL_TOTNUM][DS_SZ_KEYWORD] = { 
                            TRING_AXIS_NAM,  TRING_A_NAM,
                            TRING_B_NAM,     TRING_C_NAM, 
                            TRING_D_NAM,     TRING_BETA_NAM, 
                            TRING_GAMMA_NAM, TRING_THRESH12_NAM,
                            TRING_OFFSET_NAM } ;

         dmDescriptor *tring_cols_p[TRING_COL_TOTNUM];
         dmDataset *tap_ds_p = NULL;
         dmBlock *tap_block_p = NULL;
         char axis_val_a[DS_SZ_KEYWORD];   /* value of column AXIS */

         /* open file */
         if ((tap_ds_p = dmDatasetOpen(inp_p->tapfile)) != NULL)
               tap_block_p = dmBlockOpen(tap_ds_p, TRING_EXTNAME);

         if (tap_ds_p && tap_block_p)
         {
            /*********************************************************
             * tring_coeffs_p : a scalar that holds 2 rows of data   
             * ( tring_coeffs_p: 1 dim array with only 1 element )
             *********************************************************/
            if (((*tring_coeffs_p)=(TRING_COEFFS_P_T)calloc(1, 
                                 sizeof(TRING_COEFFS_T))) != NULL)
            {
               short col_dependency_fail = 0;

               /* open all 9 columns */
               /* add a new TRING_OFFSET column */

               set_TRING_OFFSET = 1 ;   /*get value from TRING_OFFSET column*/

               for (ii = 0; ii < TRING_COL_TOTNUM ; ii++)
               {
                  if ((tring_cols_p[ii] =
                       dmTableOpenColumn(tap_block_p,col_names_a[ii]))==NULL)
                  {
                     if (ii != (TRING_COL_TOTNUM-1))  
                        col_dependency_fail = 1;
                     else
                     {
                        /* if TRING_OFFSET column failed to open */
                        set_TRING_OFFSET = 0;   /*set TRING_OFFSET value to 0*/
                        dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                      "WARNING: %s does not contain the input column TRING_OFFSET, set value to 0.", inp_p->tapfile);
                     }
                  }
               }

               if (col_dependency_fail)  /* TRING_OFFSET not included */
               {
                   dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                      "WARNING: %s does not contain the input columns necessary to perform the tap ring test, so the test will be skipped.", inp_p->tapfile);
                   free(*tring_coeffs_p);
                   *tring_coeffs_p= NULL;
               }
               else
               {
                  /* read the two rows ; max_rows = 2 */
                  for (ii = 0; ii < max_rows; ii++)
                  {
                     memset(axis_val_a, 0, DS_SZ_KEYWORD) ;
                     dmGetScalar_c(tring_cols_p[TRING_AXIS_INDEX], axis_val_a,
                                   DS_SZ_KEYWORD-1);    /* -1 */

                     (*tring_coeffs_p)->A[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_A_INDEX]);
  
                     (*tring_coeffs_p)->B[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_B_INDEX]);

                     (*tring_coeffs_p)->C[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_C_INDEX]);

                     (*tring_coeffs_p)->D[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_D_INDEX]);

                     (*tring_coeffs_p)->BETA[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_BETA_INDEX]);

                     (*tring_coeffs_p)->GAMMA[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_GAMMA_INDEX]);

                     (*tring_coeffs_p)->THRESH12[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_THRESH12_INDEX]);

                 if (set_TRING_OFFSET == 1 )
                 {
                     (*tring_coeffs_p)->OFFSET[ii] =
                         dmGetScalar_d(tring_cols_p[TRING_OFFSET_INDEX]);
                 }
                 else
                 {
                     (*tring_coeffs_p)->OFFSET[ii] = 0.0 ;
                 }

                     /* check the value of the AXIS column : AXIS=="U" or "V"*/
                     if (strcmp(axis_val_a, TRING_AXIS_U ) == 0)
                          (*tring_coeffs_p)->U_index = ii;
                     else
                          (*tring_coeffs_p)->V_index = ii;

                     dmTableNextRow(tap_block_p);
                  } /* end: for(rows) */

               } /* end: else of if(col_dependency_fail) */

            } /* end: if(calloc) */
            else
            {
               dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Generic);
               *tring_coeffs_p= NULL;
            }
         } /* end: if (tap_ds_p && tap_block_p) */
         else
         {
            /* unable to open TRING_EXTNAME extension */
            dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                     "WARNING: Unable to open the %s extension in file %s.  The tap ring test will be skipped.",  TRING_EXTNAME, inp_p->tapfile);
         }

         if (tap_block_p) dmBlockClose(tap_block_p);
         if (tap_ds_p) dmDatasetClose(tap_ds_p);

      }
      else
      {
         dsErrAdd(hpe_err_p, dsOPENFILEWERR, Individual, Generic,
                  inp_p->tapfile);
      } /* end: if(dmDatasetAccess) */
   }
   else 
   {
      if (inp_p->debug >=1)   
         printf("Skip the tap ring test. The filename is required for this test.\n"); 
   } /* end: if(strcmp(NONE,none) etc. */
} /* end: open_tap_ring_file */


/**********************************************************************
 * check_tap_ring : 
 *
 *    In order to apply the tap ring correction, the following
 *    conditions have to be satisfied :
 *         (1)  A1, A2, A3 have to be ALL positive ;
 *         (2)  amp_sf = 3   (the corrected amp_sf from apply_amp_sf_cor)
 *         (3)  A1 > A3 
 *              ie.    for U AXIS :    AU1 > AU3 ;
 *                     for V AXIS :    AV1 > AV3 ;
 *         (4)  A1 > (THRESH12) * A2
 *              ie.    for U AXIS :    AU1 > thresh12_U * AU2 
 *                     for V AXIS :    AV1 > (thresh12_V * AV2 
 *         (5) If all TRUE, then       compute_A3_flag = True
 *              ie.    for U AXIS :    compute_AU3_flag = True 
 *                     for V AXIS :    compute_AV3_flag = True 
 *
 *    The function returns the following amplitudes :
 *         evt_p->amps_dd[plane][HDET_1ST_AMP] = raw_A1
 *         evt_p->amps_dd[plane][HDET_2ND_AMP] = raw_A2
 *         evt_p->amps_dd[plane][HDET_3RD_AMP] = tap_corrected_A3
 ***********************************************************************/
void check_tap_ring(
    INPUT_PARMS_P_T  inp_p,  /*i: obs_widthres, evt_widthres */
    EVENT_REC_P_T  evt_p,  /*u:amps_dd[][HDET_3RD_AMP],amps_tap_flag*/
    TRING_COEFFS_P_T  tring_coeffs_p  /* i */  )
{
  short  axis_idx = HDET_PLANE_X;               /*axis index; 0 for U axis */
  short  coeffs_idx = tring_coeffs_p->U_index;  /*coeff index for U axis; 
                                                  value can be 0 or 1*/
  short  amp_sf ;
  short  ii ;
  short  positive3a ;  /* 1 means (A1,A2,A3) all positive; otherwise 0 */

  short  WIDTHRES ; 

  double AMP1 = 0.0;         /* A1   */
  double AMP2 = 0.0;         /* A2   */
  double AMP3 = 0.0;         /* A3   */
  double AMP3_corr = 0.0;    /* A3'  */

  /* coefficients of one row from the tap ring file */
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double D = 0.0;
  double BETA  = 0.0;
  double GAMMA = 0.0;
  double THRESH12 = 0.0;
  double OFFSET = 0.0;    /* 4/2003 */
  double theta ;    /* sin(theta) ; in radian;*/
  double ttt = 0.0;       /* 11/2008 */

  HRC_STATUS_T  bad_bit_mask=HDET_U_TAP_RING_STS ;   /* JCC(6/16/00) */
  HRC_STATUS_T  widthExceedBit=HDET_U_WIALIGN_STS; /* long; bit 20 */

  /* get the value of amp_sf  */
  amp_sf =  evt_p->amp_sf ;  /*the corrected amp_sf from apply_amp_sf_cor*/

  /* ----------------------------------------------------------------
   * if (amp_sf==3) && (A1 > A3) && (A1 > thresh12 * A2 )
   * then  amp_tap_flag = 1    ( ie. apply the tap corrections )
   * --------------------------------------------------------------*/
   /* loop for each axis (U and V Axes), check the conditions for
    * each axis, and set the condition flags appropriate. */
   for (ii = 0; ii < 2; ii++)
   {
      /*JCC(7/2003)-reset status bit for tapRing  before reapply the corr.*/
      /* (evt_p->status&bad_bit_mask)?printf("before 1\n"):printf(""); */
       evt_p->status &= (~bad_bit_mask) ;

      /* amps_dd are RAW amps from the infile and assume they're >= 0 */
      /* modify amps_dd[u:v][HDET_3RD_AMP] for tap ring at the end*/
       AMP1 = evt_p->amps_dd[axis_idx][HDET_1ST_AMP];
       AMP2 = evt_p->amps_dd[axis_idx][HDET_2ND_AMP];
       AMP3 = evt_p->amps_dd[axis_idx][HDET_3RD_AMP];

       /* Are all of A1, A2 and A3 positive ? */
       if ((AMP1 >= 0.0)&&(AMP2 >= 0.0)&&(AMP3 >= 0.0))
          positive3a = 1 ;   /* yes, all positive */
       else 
          positive3a = 0 ;

       /* read test coefficients for U Axis and then V Axis*/
       A = tring_coeffs_p->A[coeffs_idx];
       B = tring_coeffs_p->B[coeffs_idx];
       C = tring_coeffs_p->C[coeffs_idx];
       D = tring_coeffs_p->D[coeffs_idx];
       BETA = tring_coeffs_p->BETA[coeffs_idx];
       GAMMA = tring_coeffs_p->GAMMA[coeffs_idx];
       THRESH12 = tring_coeffs_p->THRESH12[coeffs_idx];
       OFFSET = tring_coeffs_p->OFFSET[coeffs_idx];
      /*--------------------------------------------------
       *(4/2003)-take width_threshold in obs.par  or
       * WIDTHRES in evt; obs par overwrite evt1 ;
       *--------------------------------------------------*/
       WIDTHRES = NEG_9999 ;   /* initial */
       if (inp_p->obs_widthres != NEG_9999 )
           WIDTHRES =  inp_p->obs_widthres ;
       else
           WIDTHRES =  inp_p->evt_widthres ;

       /* widthExceedBit = 0 or (2**20) or (2**21) */
       widthExceedBit &= evt_p->status; 

       /*--------------------------------------------------*/
       /* check the conditions and compute the corrected A3 */
       /*--------------------------------------------------*/
       ttt = THRESH12 * AMP2 + OFFSET ;                     /* 11/2008 */
       if( (positive3a==1)&&
           (amp_sf==3) && (hpe_gt(AMP1,AMP3)==1) &&         /* Is amp1>amp3 ? */
           ( ( (WIDTHRES!=2) && (hpe_gt(AMP1,ttt)==1) )||   /* Is amp1>ttt ? */
             ( (WIDTHRES==2) && (widthExceedBit==0) )
           )
         )
       {
           evt_p->amps_tap_flag[axis_idx] = 1 ;  /*apply tap corrections*/ 

           /*JCC(6/16/00)- add status bits for the tap_ring test */
           evt_p->status |= bad_bit_mask;     /* set the bit to 1 */

           /*--------------------------------
            * compute the corrected A3.   
            *--------------------------------*/
           theta=2.0*PI_RAD*(AMP2-BETA*(pow(AMP2/AMP1,GAMMA)-1))/(C*AMP2+D) ;
           AMP3_corr = AMP3 - ((AMP2 + B )/A)* sin(theta) ;
       }
       else
       {
           evt_p->amps_tap_flag[axis_idx] = 0 ;  /* no tap corrections*/
           AMP3_corr  = AMP3 ;
       }

      /*-------------------------------------------------------------
       * The variable 'evt_p->amps_dd[u:v][HDET_3RD_AMP]' will no longer
       * mean 'raw A3'. It will represent the tap corrected A3 value, 
       * and will be used by other existing tests.
       *
       * Make sure the range of corrected A3 :  0 <= A3' <= 4095  
       *-------------------------------------------------------------*/
       evt_p->amps_dd[axis_idx][HDET_3RD_AMP] = AMP3_corr ;
       evt_p->amps_dd[axis_idx][HDET_3RD_AMP] =
               get_min(evt_p->amps_dd[axis_idx][HDET_3RD_AMP], 4095);
       evt_p->amps_dd[axis_idx][HDET_3RD_AMP] = 
               get_max(evt_p->amps_dd[axis_idx][HDET_3RD_AMP], 0);
  
      /*--------------------------------------------
       * continue checking the next axis (Y or V) 
       *-------------------------------------------*/
       axis_idx = HDET_PLANE_Y;
       coeffs_idx = tring_coeffs_p->V_index;
       bad_bit_mask = HDET_V_TAP_RING_STS ;    /* JCC(6/16/00) */
       widthExceedBit = HDET_V_WIALIGN_STS ; /* long; bit 21 */

   } /* end: for (ii) */
} /* end: check_tap_ring */
