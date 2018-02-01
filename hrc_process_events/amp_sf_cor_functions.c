/*                                                                
**  Copyright (C) 1997-2009  Smithsonian Astrophysical Observatory 
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
* FILE NAME:    amp_sf_cor_functions.c
* DEVELOPEMENT: tools
* DESCRIPTION: This file contains functions for 'amp_sf corrections'
* FUNCTIONS : 
*     open_amp_sf_cor_file() ;
*     sum_raw_amps() ;
*     apply_amp_sf_cor() ;
*     write_amp_sf_corr_key() ;
*
* Use the value of 'range_switch_level' (from obs.par or evt1) to find the matched
* row data from amp_sf_corr caldb file. Once it's found, get the rest of columns 
* 'pha_1to2,width_1to2, pha_2to3, width_2to3 and gain', and recompute 'amp_sf' 
* with the given algorithm.
*
* The modified amp_sf will be used for tap_ring correction, amp saturation
* test, pha_ratio check and will be written to the output data. 
* 
* The output header will add a new keyword AMPSFCOR
* (if 'do_amp_sf_cor==1 (yes)', AMPSFCOR = 1 ;  )
* (if 'do_amp_sf_cor==0 (no )', AMPSFCOR = 0 ;  )
*
* Notes :
*  Place the new function 'sum_raw_amps()' before tap_ring_corrections,
*  so it will return the sum('raw_sum_amps') of RAW amplitudes('amps_sh').
*
* JCC(8/02)-initial version
*          - if (evt_AMPSFCOR==false)&&(AMPSFCOR==false), then write 
*            AMPSFCOR=false to the outfile ; otherwise, set it to 'true'; 
* JCC(9/02)-use "<=" instead of "<" for diff1,2,3 in apply_amp_sf_cor(). 
*          -bug fixed in sum_raw_amps().
(2/2005)-remove unused variable
(1/2009)-updated comments.
*H***********************************************************************/
#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif

/**********************************************************************
 * open_amp_sf_cor_file :
 *   Get the column values from the caldb file of amp_sf correction; 
 **********************************************************************/
void open_amp_sf_cor_file ( 
        INPUT_PARMS_P_T    inp_p,             /*i: ampsfcorfile + debug*/
        AMPSFCOR_COEFF_P_T *ampsfcor_coeff,   /*o: col's in amp_sf_corr caldb*/
        dsErrList          *hpe_err_p         /*u */  )
{
  short ii;
  long   tmp_row=0, tot_row = 0  ; 
  long   curr_row = 0 ;
  long   row_check = dmSUCCESS;
  char  col_names_a[AMPSFCOR_COL_TOTNUM][DS_SZ_KEYWORD] = {
              AMPSFCOR_RANGE_SWITCH_LEVEL, AMPSFCOR_RANGE_SWITCH_TOL,
              AMPSFCOR_PHA_1TO2, AMPSFCOR_WIDTH_1TO2,  AMPSFCOR_PHA_2TO3,
              AMPSFCOR_WIDTH_2TO3,  AMPSFCOR_GAIN   } ;

  dmDescriptor *ampsfcor_col[AMPSFCOR_COL_TOTNUM];
  dmDataset *in_ds = NULL;
  dmBlock *inblock = NULL;

  *ampsfcor_coeff = NULL;      /* a scalr holds 2 rows of data */

  /* make sure a file name is provided before processing */
  if ((strcmp(inp_p->ampsfcorfile, "NONE") != 0) &&
      (strcmp(inp_p->ampsfcorfile, "none") != 0) &&
      (strcmp(inp_p->ampsfcorfile, "\0")   != 0))
  {
     if (dmDatasetAccess(inp_p->ampsfcorfile, "R") == dmTRUE)
     {
        /* open file */
        if ((in_ds = dmDatasetOpen(inp_p->ampsfcorfile)) != NULL)
              inblock = dmBlockOpen(in_ds, AMPSFCOR_EXTNAME);

        if (in_ds && inblock)
        {
           /*********************************************************
            * ampsfcor_coeff : 1 dimension array with 6 elements that 
            *                  holds one row of data with 6 columns. 
            *********************************************************/
           if (((*ampsfcor_coeff)=(AMPSFCOR_COEFF_P_T)calloc(1, 
                                sizeof(AMPSFCOR_COEFF_T))) != NULL)
           {
              short col_dependency_fail = 0;

             /*----------------------
              *  open all 6 columns 
              *----------------------*/
              for (ii = 0; ii < AMPSFCOR_COL_TOTNUM ; ii++)
              {
                 if ((ampsfcor_col[ii] =
                      dmTableOpenColumn(inblock,col_names_a[ii]))==NULL)
                 {
                    col_dependency_fail = 1;
                 }
              }

              if (col_dependency_fail)
              {
                 dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                    "WARNING: %s does not contain the input columns necessary to perform the amp_sf correction, so the test will be skipped.", inp_p->ampsfcorfile);
                 free(*ampsfcor_coeff);
                 *ampsfcor_coeff = NULL;
              }
              else
              {
                 /*---------------------------------------------------------
                  * loop thru # of rows in the amp_sf_corr caldb file until 
                  *  the matched RANGE_SWITCH_LEVEL is found. 
                  *--------------------------------------------------------*/
                  tmp_row = dmTableGetRowNo( inblock ) ;  /* current row numb. */
                  tot_row = dmTableGetNoRows( inblock ) ;  /* total row numb.   */
                  row_check = dmSUCCESS ;

                  if ( tmp_row != dmBADROW)
                  while ((row_check != dmNOMOREROWS) &&      /* while(evt_next_row)*/
                         (hpe_err_p->contains_fatal == 0) &&
                         (curr_row < tot_row  ) )
                  {
                     (*ampsfcor_coeff)->range_switch_level = 
                     dmGetScalar_s(ampsfcor_col[AMPSFCOR_RANGE_SWITCH_LEVEL_IDX]);

                     (*ampsfcor_coeff)->range_switch_tol = 
                     dmGetScalar_s(ampsfcor_col[AMPSFCOR_RANGE_SWITCH_TOL_IDX]);

                    /*--------------------------------------------------------------------
                     * If the difference betw the column of RANGE_SWITCH_LEVEL in caldb 
                     * and the one from obs.par (or from evt1) is smaller than the column 
                     * of RANGE_SWITCH_TOL in caldb, then get the rest columns from caldb. 
                     *
                     * Otherwise, go to the next row.
                     *-------------------------------------------------------------------*/
                     if ( abs((*ampsfcor_coeff)->range_switch_level - inp_p->range_switch_level) <= 
                           (*ampsfcor_coeff)->range_switch_tol )
                     {
                        inp_p->match_range_switch_level = TRUE ; /* found the matched range_switch_level */

                        (*ampsfcor_coeff)->PHA_1TO2 =
                              dmGetScalar_d(ampsfcor_col[AMPSFCOR_PHA_1TO2_IDX]);

                        (*ampsfcor_coeff)->WIDTH_1TO2 =
                            dmGetScalar_d(ampsfcor_col[AMPSFCOR_WIDTH_1TO2_IDX]);

                        (*ampsfcor_coeff)->PHA_2TO3 =
                            dmGetScalar_d(ampsfcor_col[AMPSFCOR_PHA_2TO3_IDX]);

                        (*ampsfcor_coeff)->WIDTH_2TO3 =
                            dmGetScalar_d(ampsfcor_col[AMPSFCOR_WIDTH_2TO3_IDX]);

                        (*ampsfcor_coeff)->GAIN =
                            dmGetScalar_d(ampsfcor_col[AMPSFCOR_GAIN_IDX]);

                        break ;    /* break the 'for' loop   */

                     }/* end: if [ abs (caldb->range_switch_level - inp_p->range_switch_level ] <=
                                 caldb->range_switch_tol) */

                     row_check = dmTableNextRow(inblock);
                     curr_row =  curr_row + 1 ;

                  } /* end :  while */
                  if ( inp_p->match_range_switch_level == FALSE )
                  {
                      dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                       "WARNING: Skip the amp_sf correction.  Can't find a matched range_switch_level in amp_sf_corr caldb file.\n") ;
                  }

               } /* end: else of if(col_dependency_fail) */
            } /* end: if(calloc) */
            else
            {
               dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Generic);
               *ampsfcor_coeff = NULL;
            }
        } /* end: if (in_ds && inblock) */
        else
        {
           /* unable to open AMPSFCOR_EXTNAME extension */
           dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
               "WARNING: Unable to open the %s extension in file %s.  The amp_sf correction will be skipped.", AMPSFCOR_EXTNAME, inp_p->ampsfcorfile);
        }

        if (inblock) dmBlockClose(inblock);
        if (in_ds) dmDatasetClose(in_ds);

     }
     else
     {
        dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
           "WARNING: Skip amp_sf correction. Unable to open its caldb file.\n") ;
     } /* end: if(dmDatasetAccess) */
  }
  else 
  {
        dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
      "WARNING: Skip the amp_sf correction. The caldb filename is required.");
  } /* end: if(strcmp(NONE,none) etc. */

  /* ---- for the output keyword : AMPSFCOR ---- */
  if (inp_p->match_range_switch_level == TRUE )
  {
      inp_p->AMPSFCOR = TRUE ;
  }
  else
  {
      inp_p->AMPSFCOR = FALSE ;
  }
} /* end: open_amp_sf_cor_file */

/**************************************************************
 * all all 6 raw amplitudes (AU1,AU2,AU3,AV1,AV2,AV3)
 * and return sum_raw_amps.
 **************************************************************/
void sum_raw_amps (
      EVENT_REC_P_T evt_p,   /* I:  amps_sh[0-1][0-2] */
      double   *tot_raw_amps /* O: */ )
{
   short   tmp = 0 ;
   short   count;
  
   *tot_raw_amps = 0.0 ;

   /*---- HDET_1ST_AMP=0; HDET_NUM_AMPS=3 ----*/
   for (count = HDET_1ST_AMP; count < HDET_NUM_AMPS; count++)
   {
      tmp  = tmp +  evt_p->amps_sh[HDET_PLANE_X][count] 
                 +  evt_p->amps_sh[HDET_PLANE_Y][count];
   }
  *tot_raw_amps = (double)tmp ;

}  /* end : sum_raw_amps */

/*****************************************************************
 * apply_amp_sf_cor() uses raw-amplitudes 'au1-3,av1-3' in infile 
 * to adjust 'evt_p->amp_sf'. The main code will then propogate it 
 * to the output column 'amp_sf'.
 *****************************************************************/
void apply_amp_sf_cor (
          EVENT_REC_P_T      evt_p,          /* u: amp_sf */
          AMPSFCOR_COEFF_P_T ampsfcor_coeff  /* i:  */  
         )
{
  short    AMP_SF = 3 ;   /* default */

  /* column data from the amp_sf_corr caldb file */
  double PHA_1TO2 = 0.0 ;            /* column PHA_1TO2 */
  double WIDTH_1TO2 = 0.0 ;          /* column WIDTH_1TO2 */
  double PHA_2TO3 = 0.0 ;            /* column PHA_2TO3 */
  double WIDTH_2TO3 = 0.0 ;          /* columnn WIDTH_2TO3 */
  double GAIN = 0.0 ;                /* column GAIN */

  double Pha ;                       /* pha column in event file */

  double   SUMAMPS = 0.0 ;    /* sum of 6 raw amplitudes */
  double  SCALED_SUMAMPS = 0.0 ; 
  double   diff1=0.0, diff2 = 0.0 , diff3 = 0.0 ;

 /*-----------------------------------------------------------------
  * Compute the sum of all 6 raw amplitude ;
  *
  * Make sure 'sum_raw_amps()' will be called event by event ;
  *
  * The returned 'SUMAMPS' is the sum of 6 raw amplitudes and will be 
  * used only for amp_sf_corr in this function;
  *----------------------------------------------------------------*/
  sum_raw_amps ( evt_p, &SUMAMPS) ;

 /* read test coefficients for U Axis and then V Axis*/
 /* RANGE_SWITCH_LEVEL = ampsfcor_coeff->RANGE_SWITCH_LEVEL ; */
  PHA_1TO2 = ampsfcor_coeff->PHA_1TO2 ;
  WIDTH_1TO2 = ampsfcor_coeff->WIDTH_1TO2 ;
  PHA_2TO3 = ampsfcor_coeff->PHA_2TO3 ;
  WIDTH_2TO3 = ampsfcor_coeff->WIDTH_2TO3 ;
  GAIN = ampsfcor_coeff->GAIN ;

 /*------------------------------------------------
  * apply the algorithm to make amp_sf corrections 
  *-----------------------------------------------*/

  Pha = ( (double)evt_p->pha ) ;

  if ( Pha < (PHA_1TO2 - WIDTH_1TO2))
  {
      AMP_SF = 1 ;
  }
  else if (Pha < (PHA_1TO2 + WIDTH_1TO2))
  {
      AMP_SF = 1 ;
      SCALED_SUMAMPS = 0.5 * SUMAMPS / GAIN ;  
      diff1 = fabs(Pha - SCALED_SUMAMPS) ;
      diff2 = fabs(Pha - 2.0 *SCALED_SUMAMPS) ;
      if (diff2 <= diff1) 
          AMP_SF = 2 ;
  }
  else if (Pha < (PHA_2TO3 - WIDTH_2TO3))
  {
      AMP_SF = 2 ;
  }
  else if (Pha < (PHA_2TO3 + WIDTH_2TO3))
  {
      AMP_SF = 2 ;
      SCALED_SUMAMPS = 0.5*SUMAMPS/GAIN ;
      diff2 = fabs(Pha - 2.0 *SCALED_SUMAMPS) ;
      diff3 = fabs(Pha - 4.0 *SCALED_SUMAMPS) ;
      if (diff3 <= diff2) 
          AMP_SF = 3 ;
  }
  else
  {
      AMP_SF = 3 ;
  }

  /* --------------------------------------------------------------
   * reassign the raw 'amp_sf' (from infile) with the corrected
   * one (ie. AMP_SF) and will be written to the outColumn 'amp_sf'.
   * ------------------------------------------------------------*/
  evt_p->amp_sf = AMP_SF ;  

} /* end: apply_amp_sf_cor */


/*H***********************************************************************
* FILE NAME: write_amp_sf_corr.c
 
* DEVELOPEMENT: tools
 
* JCC(8/2002) - write the amp_sf correction keywords to the output

   AMPSFCOR :  TRUE  = 1 = amp_sf_corr is applied ; 
               FALSE = 0 = amp_sf_corr is not applied;
   RANGELEV :  the value of range_switch_level from obs.par 
               or RANGELEV from evt1 header; default=0 
   -------
   If (inp_p->do_amp_sf_corr==TRUE)
     The output AMPSFCOR can be TRUE or FALSE ;
      ( amp_sf_cor caldb not found, output AMPSFCOR=FALSE ; )
      ( or if caldb found, but no match to 
          range_switch_level (obs.par or evt1), output  AMPSFCOR=FALSE ; )
   If (inp_p->do_amp_sf_corr==FALSE)
     the output AMPSFCOR will be FALSE
*H***********************************************************************/


#ifndef HRC_PROCESS_EVENTS_H 
#include "hrc_process_events.h"
#endif

/*************************************************************************
* DESCRIPTION
*
* For the keyword : RANGELEV
*     if get_range_switch_level==TRUE
*         write RANGELEV 
*     else  don't write it ;
* 
* For the keyword : AMPSFCOR
*     always write it ; 
*************************************************************************/
void write_amp_sf_corr_key(
   dmBlock*        extension,/*I - output file name */
   INPUT_PARMS_P_T inp_p     /*I - AMPSFCOR from the amp_sf code ;
                             range_switch_level from obs.par or evt1 header*/
   )
{
    /*-----------------------------------------------------
     * write RANGELEV only if it's found in obsfile or evt1
     * ( and regardless 'do_amp_sf_cor' )
     *-----------------------------------------------------*/
    if (inp_p->get_range_switch_level==TRUE)  /* for:  do_amp_sf_cor = yes */
    {
        dmKeyWrite_s(extension, "RANGELEV", inp_p->range_switch_level,
                      NULL, "range_switch_level") ;   /* in short dtype */
    }
    else  /*  for:  do_amp_sf_cor = no or yes */
    {
        if (inp_p->obs_r_s_l >= 0 )
        {
            dmKeyWrite_s(extension, "RANGELEV", inp_p->obs_r_s_l,
                      NULL, "range_switch_level") ;   /* in short dtype */
        }
        else
        {
            if (inp_p->evt_r_s_l >= 0 )
            {
                dmKeyWrite_s(extension, "RANGELEV", inp_p->evt_r_s_l,
                      NULL, "range_switch_level") ;   /* in short dtype */
            }
        }
    }
    /*------------------------------
     * always write out AMPSFCOR  
     *------------------------------*/
    if ((inp_p->evt_AMPSFCOR==FALSE)&&(inp_p->AMPSFCOR==FALSE))
       inp_p->AMPSFCOR = FALSE ;
    else
       inp_p->AMPSFCOR = TRUE  ;

    dmKeyWrite_q(extension, "AMPSFCOR", inp_p->AMPSFCOR,
                  NULL,"TRUE=apply_amp_sf_corr;  FALSE=NOT_apply;");

} /* end: write_amp_sf_corr_key */
