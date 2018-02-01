/*                                                                
**  Copyright (C) 1996-2007  Smithsonian Astrophysical Observatory 
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


#ifndef HOT_PIXEL_DEFS_H
#define HOT_PIXEL_DEFS_H 


#ifndef STKLIB
#include <stack.h>
#define STKLIB
#endif

#ifndef ASCDM_H
#include "ascdm.h"
#endif

#include <stdio.h>

/*  the following defines are used to map the columns of the bad pixel file
 *  to the fields of the bad pixel record structure.
 *
 *  BAD PIXEL MAPPINGS 
 */ 

#define BAD_RAWX        0
#define BAD_RAWY        1
#define BAD_CHIPID      2
#define BAD_UNKNOWN     3    


 
/*  the following defines are used for hot spot/bad pixel file columns. The 
 *  names below are hardwired based on the names that are anticipated to be 
 *  received from telemetry and may need revision. 
 *
 *  BAD PIXEL FILE COLUMNS
 */
 
#define BAD_RAWX_COL        "RAWX"
#define BAD_RAWY_COL        "RAWY"
#define BAD_CHIPID_COL      "CHIP_ID" 



/*  the following structure holds all information pertaining to any single
 *  badpixel (hot spot) record and is used to pass this data to routines 
 *  for processing.
 *
 *  BAD PIXEL RECORD STRUCTURE
 */
 
typedef struct bad_pix_t {
   int  x;    /* raw x coord */
   int  y;    /* raw y coord */
   struct bad_pix_t *next;   /* pointer to next bad pixel in list */
} BAD_PIX_T,  *BAD_PIX_P_T;

typedef BAD_PIX_P_T BAD_PIX_A_T[3];


/*  The following structure is used by hrc_process_events to store information
 *  pertaining to bad pixels. The fields in this structure are managed in
 *  functions in hot_pixel_functions.c
 *
 *  BAD PIXEL SETUP STRUCTURE
 */
typedef struct bad_pix_setup_t {
   Stack       stack;       /* input aspect file(s) 'stack'     */
   dmDataset*    dataset;  /* dataset (file) handle */
   dmBlock*      primary;  /* header keyword extension handle */
   dmBlock*      extension; /* event data extension handle */
   dmDescriptor**  desc;   /* event attributes */
   int           num_cols; /* number of columns in event file */
   int           num_rows; /* number of rows in the event file */
   int           row_check; /* indicates if last row has been read */
   int           curr_row;  /* the current row being read from file */
   dmDataType*   types;    /* ptr to array of event data column types */
   short*        mapping;  /* mapping between data structure and file columns */
   char*         file;     /* name of bad pixel file */ 
} BAD_PIX_SETUP_T, *BAD_PIX_SETUP_P_T;
 
 



#endif /* must be last line of file */
