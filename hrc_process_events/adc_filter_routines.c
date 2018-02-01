/*                                                                
**  Copyright (C) 2000-2007,2010  Smithsonian Astrophysical Observatory 
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
 
* FILE NAME: adc_filter_routines.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: This file contains several functions used by 
  hrc_process_events to evaluate an event as good or bad based on various
  adc criteria.

  The modules in this file are:

  open_amp_saturation_file()
  open_evt_flatness_file()
  open_hyperbolic_file()
  check_amp_saturation()
  check_evt_flatness()
  check_hyperbolic()  
  hyperbolic_help()

  For details on the functionality of any of the above listed modules, 
  please see the 'description' comment preceding the specific module's 
  source code.

  JCC(5/11/00)-Use "double amps_dd" to compute the corrections.
  JCC(5/12/00)-Add a new column of amp_sf to the saturation file and the 
               file will have 4 rows instead of 1. Rewrite the function
               "open_amp_saturation_file" and update the function
               "check_amp_saturation".
  JCC(5/16/00)-Add a new coeff C to the hyperbolic file. Update 
               open_hyperbolic_file() and check_hyperbolic().
  JCC(5/24/00)-accept both old and new formats for the sat and hyper files. 
               update open_amp_saturation_file() and open_hyperbolic_file().
  JCC(12/21/00)-check on nonzero division for hyperb. and flatness tests, 
                otherwise, it will give 'floating point error' on Alpha.
               -enhanced to remove warnings from compiling.
  JCC(2/9/01) - (bugfix) memset the variable before calling dmGetScalar_c
  JCC(6/7/01) - passing maxlen-1 to dmGetScalar_c
  JCC(7/2003)-reset status bit before reapply any correction, such as,
                  hyperbolic amp_saturation evt_flatness 
*H***********************************************************************/

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#define HRC_PROCESS_EVENTS_H
#endif 

/* routine to open/read hyperbolic test ARD.
 * JCC(5/16/00) - Add a new coeff C.
 * JCC(5/24/00)- the old hyper file has 2 rows and 5 columns, and 
                 the new one has 2 rows and 6 columns.
   JCC(2/9/01) - memset the variable before calling dmGetScalar_c
 */
void open_hyperbolic_file(
     char *hyp_file,
     HYP_TEST_P_T *hyp_coeffs_p,
     dsErrList *hpe_err_p)
{
   short  new_format = 1 ;   /* 1: new file ; 0: old file */
   short  max_rows = 2;      /* 2: for both new and old file  */ 

   *hyp_coeffs_p = NULL;

   /* make sure a file name is provided before processing */
   if(((strcmp(hyp_file, "NONE") != 0) && 
      (strcmp(hyp_file, "none") != 0)) &&
      (strcmp(hyp_file, "\0") != 0))
   {
      if (dmDatasetAccess(hyp_file, "R") == dmTRUE)
      {
         short ii ;
	 char col_names_a[HYPER_COL_TOTNUM][DS_SZ_KEYWORD] = {
                                               HYPERTEST_B_COL,
					       HYPERTEST_H_COL,
					       HYPERTEST_A_COL,
					       HYPERTEST_HDELTA_COL,
					       HYPERTEST_AXIS_COL,
					       HYPERTEST_C_COL };
	 dmDescriptor *hyper_cols_p[HYPER_COL_TOTNUM];
	 dmDataset *hyp_ds_p = NULL;
	 dmBlock *hyp_block_p = NULL;
	 char axis_val_a[DS_SZ_KEYWORD];

         /* open file */
	 if ((hyp_ds_p = dmDatasetOpen(hyp_file)) != NULL)
	 {
	    hyp_block_p = dmBlockOpen(hyp_ds_p, HPE_HYP_EXTNAME);
	 }

         if (hyp_ds_p && hyp_block_p)
         {
	    if (((*hyp_coeffs_p) = (HYP_TEST_P_T)calloc(1, sizeof(HYP_TEST_T)))
		!= NULL)
	    {
	       short col_dependency_fail = 0;

	       /* open the columns (0:4) */
	       for (ii=0; ii < HYPER_COL_TOTNUM - 1; ii++) /*skip HYPER_C col*/
	       {
		  if ((hyper_cols_p[ii] = 
		       dmTableOpenColumn(hyp_block_p, col_names_a[ii])) == NULL)
		  {
		     col_dependency_fail = 1;
		  }
	       }
	       /* open the the 5th column - new HYPER_C  */
               if ((hyper_cols_p[HYPER_C_INDEX] = dmTableOpenColumn(
                         hyp_block_p, col_names_a[HYPER_C_INDEX])) == NULL)
               {
                  new_format = 0 ;      /* old file; set HYPER_C=0 later */
               }
               
	       if (col_dependency_fail)     /* HYPER_C is not included */
	       {
		  dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                          "WARNING: %s does not contain the input columns necessary to perform the hyperbolic test, so the test will be skipped.", hyp_file);
		  free(*hyp_coeffs_p);
		  *hyp_coeffs_p = NULL;
	       }
	       else
	       {
		  /* read the two rows */
		  for (ii = 0; ii < max_rows; ii++)
		  {
                     memset(axis_val_a, 0, DS_SZ_KEYWORD ) ;
		     dmGetScalar_c(hyper_cols_p[HYPER_AXIS_INDEX], axis_val_a, 
				   DS_SZ_KEYWORD-1);  /* -1 */
		     (*hyp_coeffs_p)->B[ii] = 
		           dmGetScalar_d(hyper_cols_p[HYPER_B_INDEX]);
		     (*hyp_coeffs_p)->H[ii] = 
		           dmGetScalar_d(hyper_cols_p[HYPER_H_INDEX]);
		     (*hyp_coeffs_p)->A[ii] = 
		           dmGetScalar_d(hyper_cols_p[HYPER_A_INDEX]);
		     (*hyp_coeffs_p)->H_delta[ii] = 
		           dmGetScalar_d(hyper_cols_p[HYPER_HDELTA_INDEX]);

		     if (strcmp(axis_val_a, HYP_AXIS_U) == 0)
		        (*hyp_coeffs_p)->U_index = ii;
		     else
		        (*hyp_coeffs_p)->V_index = ii;

                     /* HYPER_C exists in the new sat file */
                     if (new_format==1)
                     {
		         (*hyp_coeffs_p)->C[ii] = 
		               dmGetScalar_d(hyper_cols_p[HYPER_C_INDEX]);
                     }
                     else 
                     {
                         (*hyp_coeffs_p)->C[ii] = (double)0.0 ;
                     }

		     dmTableNextRow(hyp_block_p);
		  } /* end for(rows) */

	       } /* end else of if(col_dependency_fail) */

	    } /* end if(calloc) */
	    else
	    {
	      dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Generic);
	      *hyp_coeffs_p = NULL;
	    }

         } /* end if (hyp_ds_p && hyp_block_p) */
         else
         {
           /* unable to open HPE_HYP_EXTNAME extension */ 
            dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                     "WARNING: Unable to open the %s extension in file %s.  The hyperbolic test will be skipped.",
                     HPE_HYP_EXTNAME, hyp_file);
         }
	 
	 if (hyp_block_p) dmBlockClose(hyp_block_p);
	 if (hyp_ds_p) dmDatasetClose(hyp_ds_p);

      } /* end if(dmDatasetAccess) */
      else
      {
	 dsErrAdd(hpe_err_p, dsOPENFILEWERR, Individual, Generic,
		  hyp_file); 
      }

   } /* end if(strcmp(NONE) etc. */
} /*end: open_hyperbolic_file */

/* routine to open/read flatness test ARD */
void open_evt_flatness_file(
     char *flatness_file,
     double **flatness_coeff_p,
     dsErrList *hpe_err_p)
{
   *flatness_coeff_p = NULL;

   /* make sure a file name is provided before processing */
   if(((strcmp(flatness_file, "NONE") != 0) && 
      (strcmp(flatness_file, "none") != 0)) &&
      (strcmp(flatness_file, "\0") != 0))
   {
      if (dmDatasetAccess(flatness_file, "R") == dmTRUE)
      {
	 dmDescriptor *flatness_col_p;
	 dmDataset *flat_ds_p = NULL;
	 dmBlock *flat_block_p = NULL;

         /* open file */
         if ((flat_ds_p = dmDatasetOpen(flatness_file)) != NULL)
         {
	    flat_block_p = dmBlockOpen(flat_ds_p, HPE_AMP_FLAT_EXTNAME);
	 }

	 if (flat_ds_p && flat_block_p)
	 {
	    if (((*flatness_coeff_p) = (double *)calloc(1, sizeof(double)))
		!= NULL)
	    {
	       /* open the column */
	       if ((flatness_col_p = dmTableOpenColumn(flat_block_p, 
						       FLATTEST_LIMIT_COL))
		   != NULL)
	       {
		  *(*flatness_coeff_p) = dmGetScalar_d(flatness_col_p);
	       }
	       else
	       {
		  dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
			   "WARNING: %s does not contain the input columns necessary to perform the event flatness test, so the test will be skipped.", flatness_file);
		  free(*flatness_coeff_p);
		  *flatness_coeff_p = NULL;
	       }
	    } /* end if(calloc) */
	    else
	    {
	       dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Generic);
	       *flatness_coeff_p = NULL;
	    }

	 } /* end if (flat_ds_p && flat_block_p) */
	 else
	 {
	   /* unable to open HPE_AMP_FLAT_EXTNAME extension */ 
            dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                     "WARNING: Unable to open the %s extension in file %s.  The event flatness test will be skipped.",
                     HPE_AMP_FLAT_EXTNAME, flatness_file);
	 }

	 if (flat_block_p) dmBlockClose(flat_block_p);
	 if (flat_ds_p) dmDatasetClose(flat_ds_p);

      } /* end if(dmDatasetAccess) */
      else
      {
	 dsErrAdd(hpe_err_p, dsOPENFILEWERR, Individual, Generic,
		  flatness_file); 
      }

   } /* end if(strcmp(NONE) etc. */
}

/********************************************************************
 * routine to open/read saturation test ARD  
 * JCC(5/12/00) - Rewrite it. Sat coeffs are a function of amp_sf. 
 * JCC(5/24/00) - accept both old and new formats.
 *********************************************************************/
void open_amp_saturation_file(
          INPUT_PARMS_P_T inp_p,        /*i: passing ampsatfile & debug */
          SAT_TEST_P_T    *sat_coeffs_p, /*o: saturation coeffs */
          dsErrList       *hpe_err_p    /*u*/ )
{

   short amp_sf[4]={0,1,2,3} ;   /* new AMP_SF column */
   short new_format = 1 ;        /* 1=new format; 0=old format */
   short max_rows = 4;           /* 4=new format; 1=old format */

   /* make sure a file name is provided before processing */
   if(((strcmp(inp_p->ampsatfile, "NONE") != 0) &&
      (strcmp(inp_p->ampsatfile, "none") != 0)) &&
      (strcmp(inp_p->ampsatfile, "\0") != 0))
   {
      if (dmDatasetAccess(inp_p->ampsatfile, "R") == dmTRUE)
      {
         short ii, jj ;
         char  col_names_a[SAT_COL_TOTNUM][DS_SZ_KEYWORD] = {
                          SATTEST_AMPSP_NAM, SATTEST_NTAPS_NAM,
                          SATTEST_LOW_NAM,   SATTEST_HIGH_NAM } ;

         dmDescriptor *sat_cols_p[SAT_COL_TOTNUM];
         dmDataset *sat_ds_p = NULL;
         dmBlock *sat_block_p = NULL;

         /* open file */
         if ((sat_ds_p = dmDatasetOpen(inp_p->ampsatfile)) != NULL)
               sat_block_p = dmBlockOpen(sat_ds_p, HPE_AMP_SAT_EXTNAME);

         if (sat_ds_p && sat_block_p)
         {
            /*********************************************************
             * sat_coeffs_p : a scalar that holds 4 rows of data
             * ( sat_coeffs_p : 1 dim array with only 1 element )
             *********************************************************/
            if (((*sat_coeffs_p)=(SAT_TEST_P_T)calloc(1,
                                 sizeof(SAT_TEST_T))) != NULL)
            {
               short col_dependency_fail = 0;

               /*---------------------------------------------------------
                * open all columns
                *JCC(5/24/00): accept both new and old formats.
                * - The new sat file has a new column 'AMP_SF' and 4 rows.
                *---------------------------------------------------------*/ 
               if ((sat_cols_p[0] =           /* only for AMP_SF column */
                       dmTableOpenColumn(sat_block_p,col_names_a[0]))==NULL)
               {
                   new_format = 0 ;    /*old sat file has 1 row & 3 cols*/
                   max_rows  = 1 ;
                   if (inp_p->debug > 1)
                      printf("sat test file is in old format\n") ;
               }

               for (ii=1; ii<SAT_COL_TOTNUM; ii++) /*skip the AMP_SF col*/
               {
                  if ((sat_cols_p[ii] =
                       dmTableOpenColumn(sat_block_p,col_names_a[ii]))==NULL)
                  {
                     col_dependency_fail = 1;
                  }
               }
               if (col_dependency_fail)  /* AMP_SF col is not included */
               {
                  dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                  "WARNING: %s does not contain the input columns necessary to perform the saturation test, so the test will be skipped.", inp_p->ampsatfile);
                  free(*sat_coeffs_p);
                  *sat_coeffs_p = NULL;
               }
               else
               {
                  /*-------------------------------------------------*/
                  /* max_rows can be 4 for a new file, or 1 for old. */
                  /*-------------------------------------------------*/
                  for (ii = 0; ii < max_rows; ii++)
                  {
                    /*---------------------------------------------------
                     * sat_coeffs_p->ampsf_index[2]=0 : means the value 
                     * of amp_sf in the 0th row is 2 
                     *---------------------------------------------------*/
                     if (new_format==1)     /*else, use the default of amp_sf*/
                     {
                        amp_sf[ii]=dmGetScalar_s(sat_cols_p[SAT_AMPSF_INDEX]);
                        (*sat_coeffs_p)->ampsf_index[amp_sf[ii]] = ii ;
                     }

                     (*sat_coeffs_p)->ntaps[ii] =
                         dmGetScalar_s(sat_cols_p[SAT_NTAPS_INDEX]);

                     (*sat_coeffs_p)->adc_low[ii] =
                         dmGetScalar_s(sat_cols_p[SAT_LOW_INDEX]);

                     (*sat_coeffs_p)->adc_high[ii] =
                         dmGetScalar_s(sat_cols_p[SAT_HIGH_INDEX]);

                     dmTableNextRow(sat_block_p);
                  } /* end: for (ii=0; ii< max_rows) */

                  /*-----------------------------------------------------
                   * if the file is in old format, copy the 0th row to 
                   * the other 3 and set ampsf_index[0:3]={0,1,2,3}.
                   *-----------------------------------------------------*/
                  if (new_format==0)
                  {
                     (*sat_coeffs_p)->ampsf_index[0] = 0 ;

                     for (jj= 1; jj < 4 ; jj++ )
                     {
                         (*sat_coeffs_p)->ampsf_index[jj] = jj ;
 
                         (*sat_coeffs_p)->ntaps[jj] =
                                        (*sat_coeffs_p)->ntaps[0] ;
                         (*sat_coeffs_p)->adc_low[jj] =
                                        (*sat_coeffs_p)->adc_low[0] ;
                         (*sat_coeffs_p)->adc_high[jj] =
                                        (*sat_coeffs_p)->adc_high[0] ;
                     }
                  }   /* end: if (new_format==0) */

               }  /* end: else of if(col_dependency_fail) */
            }  /* end: if(calloc) */
            else
            {
               dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Generic);
               *sat_coeffs_p = NULL;
            }
         }  /* end: if (sat_ds_p && sat_block_p) */
         else
         {
            /* unable to open HPE_AMP_SAT_EXTNAME extension */
            dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                "WARNING: Unable to open the %s extension in file %s.  The saturation test will be skipped.",  HPE_AMP_SAT_EXTNAME, inp_p->ampsatfile);
         }

         if (sat_block_p) dmBlockClose(sat_block_p);
         if (sat_ds_p) dmDatasetClose(sat_ds_p);
      }
      else
      {
         dsErrAdd(hpe_err_p, dsOPENFILEWERR, Individual, Generic,
                  inp_p->ampsatfile);
      } /* end: if(dmDatasetAccess) */
   }
   else
   {
      if (inp_p->debug >=1)
        printf("Skip the saturation test. Require a filename for this test.\n");
   } /* end: if(strcmp(NONE,none) etc. */
} /* end: open_amp_saturation_file */


/* routine to evaluate hyperbolic test */
void check_hyperbolic(
     EVENT_REC_P_T evt_p,
     HYP_TEST_P_T hyp_coeffs_p)
{
  short ctr = 0;
  double H_plus_delta = 0.0;
  double H_minus_delta = 0.0;
  double H_dminus_A = 0.0;
  double H_dplus_A = 0.0;
  double sum_ADC = 0.0;
  double A = 0.0;
  double B = 0.0;
  double H = 0.0;
  double Hdelta = 0.0;
  double C = 0.0;
  double fp = 0.0;
  double fb = 0.0;
  double ADC1 = 0.0;
  double ADC2 = 0.0;
  double ADC3 = 0.0;
  HRC_STATUS_T bad_bit_mask = HDET_U_HYP_STS;
  short axis_index = HDET_PLANE_X;
  short coeffs_index = hyp_coeffs_p->U_index;

  /* loop for each axis : beginning with X/U axis */
  for (ctr = 0; ctr < 2; ctr++)
  {
     /*JCC(7/2003)-reset status bit for hyperbolic before reapply the corr.*/
     evt_p->status &= (~bad_bit_mask) ;

     /* read amp values */
     ADC1 = evt_p->amps_dd[axis_index][HDET_1ST_AMP];
     ADC2 = evt_p->amps_dd[axis_index][HDET_2ND_AMP];
     ADC3 = evt_p->amps_dd[axis_index][HDET_3RD_AMP];

     sum_ADC = ADC1 + ADC2 + ADC3;
     
     /* read test coefficients */
     A = hyp_coeffs_p->A[coeffs_index];
     B = hyp_coeffs_p->B[coeffs_index];
     H = hyp_coeffs_p->H[coeffs_index];
     Hdelta = hyp_coeffs_p->H_delta[coeffs_index];
     C = hyp_coeffs_p->C[coeffs_index];

     H_plus_delta = H + Hdelta;
     H_minus_delta = H - Hdelta;
     H_dminus_A = H_minus_delta - A;
     H_dplus_A = H_plus_delta - A;

/*----------------------------------------------------------------------
 * JCC(12/21/00) - check the division value (sum_ADC != 0.0 ) to avoid
 *                 overflow on Alpha. 
 *       ie. add :   'if ( fabs(sum_ADC) > DBL_EPSILON )'
 *----------------------------------------------------------------------*/
if ( fabs(sum_ADC) > DBL_EPSILON )      /* division 'sum_ADC' can't be zero */
{
     fp = fabs((ADC3 - ADC1) / sum_ADC);
     fb = ADC2 / sum_ADC;

       /* JCC(5/16/00) - add the new coeff C */
       /*              - fb >=0 , coeff_C >=0, fp >=0 */
       
       if (( fb < C ) || ( fb < 0.0 ))
         {
           evt_p->status |= bad_bit_mask;
         }
       else if ( fb > H_dplus_A )
	 {
	   evt_p->status |= bad_bit_mask;
	 }
       /*---------------------------------------------------------------
        * JCC(12/21/00)-add the condition 'A != 0.0' to avoid overflow 
        *            on Alpha (ie.  A is 'division' in hyperbolic_help) 
        * 
        *  ie. insert: ( fabs(A) >= DBL_EPSILON )
        *---------------------------------------------------------------*/
       else if ( (fb <= H_dminus_A ) && ( fabs(A) >= DBL_EPSILON ) && !         
		 ((fp > hyperbolic_help(fb, H_minus_delta, A, B)) && 
		  (fp < hyperbolic_help(fb, H_plus_delta,  A, B)) ) )
	 {
	   evt_p->status |= bad_bit_mask;
	 }
       /*---------------------------------------------------------------
        * JCC(12/21/00)-add the condition 'A != 0.0' to avoid overflow 
        *            on Alpha (ie.  A is 'division' in hyperbolic_help) 
        * 
        *  insert:  ( fabs(A) >= DBL_EPSILON )
        *--------------------------------------------------------------*/
       else if ( (fb > H_dminus_A ) && ( fabs(A) >= DBL_EPSILON ) && 
		  ( fp >= hyperbolic_help(fb, H_plus_delta,  A, B))  )
	 {
	   evt_p->status |= bad_bit_mask;
	 }
       
} /*--------- end:  if ( fabs(sum_ADC) > DBL_EPSILON ) ------------*/ 

     /* go to the next Y/V axis  */
     axis_index = HDET_PLANE_Y;
     bad_bit_mask = HDET_V_HYP_STS;
     coeffs_index = hyp_coeffs_p->V_index;

  } /* end: loop U,V axes:   for (ctr = 0; ctr < 2; ctr++)  */
} /* end: check_hyperbolic */



/* helper function for the hyperbolic test */
double hyperbolic_help(
       double x, 
       double y, 
       double A, 
       double B)
{
  return (B * (sqrt(pow(((x - y)/A),2.0) - 1)));
}
  

/*---------------------------------------------------------
 * routine to evaluate saturation test.
 * JCC(5/12/00)- change sat_test_coeffs_p to a structure ; 
 *             - sat coeffs are a function of amp_sf.
 *---------------------------------------------------------*/
void check_amp_saturation(
     EVENT_REC_P_T evt_p,            /* i: amp_sf */
     SAT_TEST_P_T  sat_test_coeffs_p  /*i */ )
{
  short ctr = 0;
  short axis_index = HDET_PLANE_X;
  double ADC[3];
  HRC_STATUS_T bad_bit_mask = HDET_U_AMP_SAT_STS;  
  short num_amps=0 ; 
  short amp_sf=0, coeff_index=0 ;

  amp_sf = evt_p->amp_sf ; 
  coeff_index = sat_test_coeffs_p->ampsf_index[amp_sf] ; 
  num_amps = sat_test_coeffs_p->ntaps[coeff_index] ; 

  /* loop for each axis: U and V  */
  for (ctr = 0; ctr < 2; ctr++)
  {
     short num_above = 0;
     short num_below = 0;
     short ctr2 = 0;

     /*JCC(7/2003)-reset status bit for amp_saturat. before reapply the corr*/
     evt_p->status &= (~bad_bit_mask) ;

     /* read amp values */
     ADC[0] = evt_p->amps_dd[axis_index][HDET_1ST_AMP];
     ADC[1] = evt_p->amps_dd[axis_index][HDET_2ND_AMP];
     ADC[2] = evt_p->amps_dd[axis_index][HDET_3RD_AMP];

     /*-----------------------------------------------------
      * OLD sat columns: ntaps, adc_low, adc_high; ONE ROW
      * NEW sat columns: amp_sf, ntaps, adc_low, adc_high; FOUR ROWS
      *
      * for sat_coeffs_p->ampsf_index[2] = 0 , it means the 
      * value of amp_sf in 0th row is 2 
      *-----------------------------------------------------*/
     for (ctr2 = 0; ctr2 < 3; ctr2++)
     {

       if (ADC[ctr2] < sat_test_coeffs_p->adc_low[coeff_index])
       {
	 num_below++;
       }
       else if (ADC[ctr2] > sat_test_coeffs_p->adc_high[coeff_index])
       {
	 num_above++;
       }
     }
    
     /* perform check */
     if (num_above >= num_amps || num_below >= num_amps)
     {
       evt_p->status |= bad_bit_mask;
     }
     
     axis_index = HDET_PLANE_Y;
     bad_bit_mask = HDET_V_AMP_SAT_STS;
  }    

}

/* routine to evaluate flatness test */
void check_evt_flatness(
     EVENT_REC_P_T evt_p,
     double flat_coeff)
{
  short ctr = 0;
  short axis_index = HDET_PLANE_X;
  double ADC1 = 0.0;
  double ADC2 = 0.0;
  double ADC3 = 0.0;
  HRC_STATUS_T bad_bit_mask = HDET_U_AMP_FLAT_STS;

  /* loop for each axis */
  for (ctr = 0; ctr < 2; ctr++)
  {
     /*JCC(7/2003)-reset status bit for evt_flat before reapply the corr.*/
     evt_p->status &= (~bad_bit_mask) ;

     /* read amp values */
     ADC1 = evt_p->amps_dd[axis_index][HDET_1ST_AMP];
     ADC2 = evt_p->amps_dd[axis_index][HDET_2ND_AMP];
     ADC3 = evt_p->amps_dd[axis_index][HDET_3RD_AMP];

     /* perform check */
     /* JCC(12/21/00) - ADC2 can't be zero on Alpha */
     if ( fabs(ADC2) > DBL_EPSILON)
     { 
        if ( ((ADC1/ADC2) >= flat_coeff) &&
	     ((ADC3/ADC2) >= flat_coeff) )
        {
          evt_p->status |= bad_bit_mask;
        }
     }

     axis_index = HDET_PLANE_Y;
     bad_bit_mask = HDET_V_AMP_FLAT_STS;
  }    
}
