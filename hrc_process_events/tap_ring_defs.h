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
 * JCC (5/1/00) - initial version
 *
 * This file defines the column information from the tap ring fits file.
 *
 * The tap ring file has 2 rows;  The 1st row (AXIS='U') contains all 
 * coeff's for U-axis ; The 2nd row (AXIS='V') has all coeff's for V-axis. 
 *
 *JCC(4/2003)-new tapRing spec ; add a column TRING_OFFSET ;
 ***************************************************************************/
#ifndef TAP_RING_DEFS_H
#define TAP_RING_DEFS_H

/* EXTNAME in the tap ring file */
#define  TRING_EXTNAME      "AXAF_TAPRINGTEST"

/* COLUMN NAMES in the tap ring file that defines the coefficients
 * of the tap ring test.  */

#define TRING_AXIS_NAM      "AXIS"
#define TRING_A_NAM         "TRING_A"
#define TRING_B_NAM         "TRING_B"
#define TRING_C_NAM         "TRING_C"
#define TRING_D_NAM         "TRING_D"
#define TRING_BETA_NAM      "TRING_BETA"
#define TRING_GAMMA_NAM     "TRING_GAMMA"
#define TRING_THRESH12_NAM  "TRING_THRESH12"
#define TRING_OFFSET_NAM    "TRING_OFFSET"   

/* COLUMN NUMBERS  */
#define TRING_COL_TOTNUM         9       /* tot. col. num. of tap ring file*/
#define TRING_AXIS_INDEX         0       /* 1st column */
#define TRING_A_INDEX            1       /* 2nd column */
#define TRING_B_INDEX            2
#define TRING_C_INDEX            3
#define TRING_D_INDEX            4
#define TRING_BETA_INDEX         5
#define TRING_GAMMA_INDEX        6
#define TRING_THRESH12_INDEX     7
#define TRING_OFFSET_INDEX       8       

/* the values of the AXIS column  */
#define TRING_AXIS_U      "U"
#define TRING_AXIS_V      "V"

/* COLUMN VALUES : */
/* structure that holds the coefficients for the tap ring test.  */
/* there should be only 2 rows in the fits file (U & V).         */
typedef struct tring_coeffs_t {
  double A[2];           /* A coeff values for U and V */
  double B[2];           /* B coeff values for U and V */
  double C[2];           /* C coeff values for U and V */
  double D[2];           /* D coeff values for U and V */
  double BETA[2];        /* BETA coeff values for U and V */
  double GAMMA[2];       /* GAMMA coeff values for U and V */
  double THRESH12[2];    /* (AU1/AU2) and (AV1/AV2) for U and V */
  double OFFSET[2];      /* (4/2003) */
  short U_index;         /* U_index=0 means 'the 1st-row' of AXIS is "U" */
  short V_index;         /* V_index=1 means 'the 2nd-row' of AXIS is "V" */
} TRING_COEFFS_T, *TRING_COEFFS_P_T ; 

/*-------------------------------*/
/* define some simple math stuff */
/*-------------------------------*/
#define  PI_RAD    3.1415926535897932385  /*pi in radians; ds_constant.c */
#define  get_min(a,b)  ((a) <= (b) ? (a) : (b))
#define  get_max(a,b)  ((a) >= (b) ? (a) : (b))

#endif
