/*                                                                
**  Copyright (C) 2000-2007  Smithsonian Astrophysical Observatory 
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

/************************************************************************/
/*
 * JCC(5/12/00) - add a new column of amp_sf to the saturation file
 *                and the file will have 4 rows instead of 1.
 * JCC(5/16/00) - add a new coeff C to the hyperbolic file.
 * JCC(5/24/00) - accept both old and new formats for the sat and hyper files.
 ************************************************************************/
#ifndef ADCFILTER_DEFS_H
#define ADCFILTER_DEFS_H

/* the following is needed for the hyperbolic test */
#include <math.h>

/* ADC filtering ARD column definitions. The following are the column   
 * names in the ARD files that define the coefficients of the various
 * ADC filtering tests: the hyperbolic test, the amp saturation test,
 * and the event flatness test.
 */

#define HYPERTEST_B_COL		"HYPER_B"
#define HYPERTEST_H_COL		"HYPER_H"
#define HYPERTEST_A_COL		"HYPER_A"
#define HYPERTEST_HDELTA_COL	"HYPER_DELTAH"
#define HYPERTEST_AXIS_COL	"AXIS"
#define HYPERTEST_C_COL	        "HYPER_C"    /* new column */

#define SATTEST_AMPSP_NAM       "AMP_SF"     /* new column */
#define SATTEST_NTAPS_NAM	"NTAPS"
#define SATTEST_LOW_NAM		"ADC_LOWER_LIMIT"
#define SATTEST_HIGH_NAM	"ADC_UPPER_LIMIT"

#define FLATTEST_LIMIT_COL	"FLATNESS_LIMIT"

/* the following defines are for lookup of hyperbolic test columns */
#define HYPER_B_INDEX		0
#define HYPER_H_INDEX		1
#define HYPER_A_INDEX		2
#define HYPER_HDELTA_INDEX	3
#define HYPER_AXIS_INDEX	4	
#define HYPER_C_INDEX	        5    /* new column */
#define HYPER_COL_TOTNUM        6

/* miscellaneous defines */
#define HYP_AXIS_U	"U"
#define HYP_AXIS_V	"V"

#define SAT_AMPSF_INDEX         0     /* new column */
#define SAT_NTAPS_INDEX         1
#define SAT_LOW_INDEX           2
#define SAT_HIGH_INDEX          3 
#define SAT_COL_TOTNUM          4

/******************************************************
 * the coefficients structure for the hyperbolic test
 *****************************************************/
typedef struct hyp_test_t {
  double B[2];  /* B coeff values for U and V */
  double H[2];  /* H coeff values for U and V */
  double A[2];	/* A coeff values for U and V */
  double H_delta[2]; /* H delta values for U and V */
  double C[2];  /* C coeff values for U and V */
  short U_index; /* index for U axis */
  short V_index; /* index for V axis */
} HYP_TEST_T, *HYP_TEST_P_T;

/******************************************************
 * the coefficients structure for the saturation test
 ******************************************************/
typedef struct sat_test_t {
  short ampsf_index[4]; /* ampsf_index[2]=0 :value of amp_sf in 0th row is 2*/
  short ntaps[4];       /* NTAPS coeff values */
  short adc_low[4];     /* ADC_LOWER_LIMIT coeff values */
  short adc_high[4];    /* ADC_UPPER_LIMIT coeff values */
} SAT_TEST_T, *SAT_TEST_P_T;

/** function prototypes	- these functions open and read the ARDs that
  contain the coefficients for the ADC filtering tests, and execute
  those tests. They are defined in hrc_process_events.h.  
extern void open_amp_saturation_file() ;
extern void open_evt_flatness_file(char*, double**, dsErrList*);
extern void open_hyperbolic_file(char*, HYP_TEST_P_T*, dsErrList*);
extern void check_amp_saturation(EVENT_REC_P_T, SAT_TEST_P_T); 
extern void check_evt_flatness(EVENT_REC_P_T, double);
extern void check_hyperbolic(EVENT_REC_P_T, HYP_TEST_P_T);
extern double hyperbolic_help(double, double, double, double);
************/
#endif
