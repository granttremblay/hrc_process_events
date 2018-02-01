/*                                                                
**  Copyright (C) 1997-2007  Smithsonian Astrophysical Observatory 
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

/***************************************************************************
 * JCC (8/2002) - initial version
 *
 * This file defines the column information from the amp_sf_corr fits file.
 *
 * The amp_sf_corr caldb file has the following columns :
 *     RANGE_SWITCH_LEVEL
 *     PHA_1TO2
 *     WIDTH_1TO2
 *     PHA_2TO3
 *     WIDTH_2TO3
 *     GAIN
 ***************************************************************************/
#ifndef AMP_SF_COR_H
#define AMP_SF_COR_H

/* EXTNAME in the amp_sf_cor file */
#define  AMPSFCOR_EXTNAME      "AXAF_AMP_SF_COR"

/* COLUMN NAMES in the amp_sf_corr file that defines the coefficients
 * of the amp_sf_corr test.  */

#define AMPSFCOR_RANGE_SWITCH_LEVEL   "RANGE_SWITCH_LEVEL"
#define AMPSFCOR_RANGE_SWITCH_TOL     "RANGE_SWITCH_TOL"
#define AMPSFCOR_PHA_1TO2             "PHA_1TO2"
#define AMPSFCOR_WIDTH_1TO2           "WIDTH_1TO2"
#define AMPSFCOR_PHA_2TO3             "PHA_2TO3"
#define AMPSFCOR_WIDTH_2TO3           "WIDTH_2TO3"
#define AMPSFCOR_GAIN                 "GAIN"

/* COLUMN NUMBERS  */
#define AMPSFCOR_COL_TOTNUM       7       /* tot col num of amp_sf_corr file*/
#define AMPSFCOR_RANGE_SWITCH_LEVEL_IDX 0  /* 1st column */
#define AMPSFCOR_RANGE_SWITCH_TOL_IDX 1    /* 2nd column */
#define AMPSFCOR_PHA_1TO2_IDX     2       /* 2nd column */
#define AMPSFCOR_WIDTH_1TO2_IDX   3       /* 3rd column */
#define AMPSFCOR_PHA_2TO3_IDX     4
#define AMPSFCOR_WIDTH_2TO3_IDX   5
#define AMPSFCOR_GAIN_IDX         6

/* COLUMN VALUES in the amp_sf_corr CALDB file */
typedef struct ampsfcor_coeff_str {
  short  range_switch_level ; /* column RANGE_SWITCH_LEVEL */
  short  range_switch_tol   ; /* column RANGE_SWITCH_TOL */
  double PHA_1TO2;            /* column PHA_1TO2 */ 
  double WIDTH_1TO2;          /* column WIDTH_1TO2 */
  double PHA_2TO3;            /* column PHA_2TO3 */ 
  double WIDTH_2TO3;          /* columnn WIDTH_2TO3 */
  double GAIN;                /* column GAIN */
} AMPSFCOR_COEFF_T, *AMPSFCOR_COEFF_P_T; 

#endif
