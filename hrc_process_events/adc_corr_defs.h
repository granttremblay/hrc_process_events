/*                                                                
**  Copyright (C) 1998-2007  Smithsonian Astrophysical Observatory 
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


#ifndef ADC_CORR_DEFS_H
#define ADC_CORR_DEFS_H 



/*  the following defines are used by hrc_process_events to specify various
 *  values used by the adc correction tables such as the names of the columns
 *  in the adc file.
 *
 *  ADC TABLE DEFINES
 */

#define ADC_UNKNOWN_COL      0        /* number for unknown adc column   */
#define ADC_TAP_COL          1        /* tap number                      */
#define ADC_P1_COL           2        /* ADC1 intercept value            */
#define ADC_Q1_COL           3        /* ADC1 slope value                */
#define ADC_P2_COL           4        /* ADC2 intercept value            */
#define ADC_Q2_COL           5        /* ADC2 slope value                */
#define ADC_P3_COL           6        /* ADC3 intercept value            */
#define ADC_Q3_COL           7        /* ADC3 slope value                */
#define ADC_AXIS_COL         8        /* specifies x or y axis           */
#define ADC_NUM_COLS         9        /* total num columns in adc table  */

#define ADC_AXIS_NAM         "AXIS"   /* name of axis col in adc table   */
#define ADC_TAP_NAM          "TAPNUM" /* name of axis col in adc table   */
#define ADC_P1_NAM           "P1"     /* intercept col name in adc table */
#define ADC_Q1_NAM           "Q1"     /* slope col name in adc table     */
#define ADC_P2_NAM           "P2"     /* intercept col name in adc table */
#define ADC_Q2_NAM           "Q2"     /* slope col name in adc table     */
#define ADC_P3_NAM           "P3"     /* intercept col name in adc table */
#define ADC_Q3_NAM           "Q3"     /* slope col name in adc table     */



/*  The following structure is used by hrc_process_events to store information
 *  pertaining to the adc file. The fields in this structure are managed in
 *  functions in adc_corr_routines.c
 *
 *  ADC FILE SETUP STRUCTURE
 */

typedef struct adc_setup_t {
   dmDataset*    dataset;      /* dataset (file) handle                */
   dmBlock*      primary;      /* header keyword extension handle      */
   dmBlock*      extension;    /* ADCCORR data extension handle        */
   dmDescriptor**  desc;       /* file data attributes                 */
   int           num_cols;     /* number of columns in the adc file    */
   int           num_rows;     /* number of rows in the adc file       */
   int           row_check;    /* indicates if last row has been read  */
   int           curr_row;     /* the current row being read from file */
   dmDataType*   types;        /* ptr to array of data column types    */
   short*        mapping;      /* maps data structure and file columns */
} ADC_SETUP_T, *ADC_SETUP_P_T;
 

 
/*
 *  the following typedef defines the structure used to house
 *  adc correction factors which may be loaded in from a user 
 *  provided adc file.
 *
 *  ADC TABLE ELEMENT STRUCTURE
 */

typedef struct adc_coor_t {
    float  p1;          /* intercept for ADC1       */ 
    float  q1;          /* slope for ADC1           */ 
    float  p2;          /* intercept for ADC2       */ 
    float  q2;          /* slope for ADC2           */ 
    float  p3;          /* intercept for ADC3       */ 
    float  q3;          /* slope for ADC3           */ 
    short  tap_num;     /* corresponding tab number */
} ADC_CORR_T, *ADC_CORR_P_T;
 


/*
 *  the following externs are function prototypes of the adc correction
 *  routines which have public access from other routines.
 *
 *  FUNCTION PROTOTYPES
 */

/* routine to allocate memory for adc tables */
extern boolean allocate_adc_table(INPUT_PARMS_P_T, 
                                  ADC_CORR_P_T*, 
                                  ADC_CORR_P_T*,
                                  dsErrList*);
 
/* routine to load values into adc table */
extern void adc_table_load(INPUT_PARMS_P_T, 
                              ADC_CORR_P_T, 
                              ADC_CORR_P_T,
                              dsErrList*);
 
/* routine to free dynamic memory used by adc tables */
extern void deallocate_adc_table(ADC_CORR_P_T*, 
                                 ADC_CORR_P_T*);
 
/* function to perform adc corrections */
extern void apply_adc_correction(ADC_CORR_P_T, 
                                 ADC_CORR_P_T, 
                                 INPUT_PARMS_P_T, 
                                 EVENT_REC_P_T);

#endif   /* last line of header file- closes #ifndef ADC_CORR_DEFS_H */
