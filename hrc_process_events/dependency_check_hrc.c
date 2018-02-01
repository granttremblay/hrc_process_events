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

/*H***********************************************************************
 
* FILE NAME: dependency_check_hrc.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: The routine dependency_check_hrc is called upon by 
  hrc_process_events to ensure that all event columns that are needed for 
  detector coordinate calculations exist in the input file. The routine 
  does not validate the data but merely verifies that the appropriate 
  columns exist. The routine dependency_check_init is called to set up a 
  bit mask used to determine if all required data dependencies have been 
  fullfilled. 
 
* NOTES:

  wmclaugh@cfa	Apr 04, 1996  First Version.

* REVISION HISTORY:
 
        Ref. No.        Date
        --------        ----
        1.1             04 Apr 1996
*H***********************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

/*H***********************************************************************
 
* DESCRIPTION: The routine dependency_check_init is called on by 
  hrc_process_events to set up the data dependency bit mask that will be
  used by dependency_check_hrc. The routine sets an appropriate bit for  
  each existing column in the input data which may be part of a necessary 
  data dependency.       

**************************************************************************/

unsigned short dependency_check_init(
   short* evt_map,      /* I - maps input file cols to event record      */  
   int    columns )     /* I - number of columns in event                */
{
   int count; 
   unsigned short mask = HDET_MASK_INIT; 

   for (count = 0; count < columns; count++)
   {
      switch (evt_map[count])
      {
         case HDET_CP_X: 
            mask |= HDET_MASK_CPX; 
         break; 

         case HDET_CP_Y: 
            mask |= HDET_MASK_CPY; 
         break; 

         case HDET_AX_1: 
            mask |= HDET_MASK_AX1; 
         break; 

         case HDET_AX_2: 
            mask |= HDET_MASK_AX2; 
         break; 

         case HDET_AX_3: 
            mask |= HDET_MASK_AX3; 
         break; 

         case HDET_AY_1: 
            mask |= HDET_MASK_AY1; 
         break; 

         case HDET_AY_2: 
            mask |= HDET_MASK_AY2; 
         break; 

         case HDET_AY_3: 
            mask |= HDET_MASK_AY3; 
         break; 

         case HDET_CHIP_X:
            mask |= HDET_MASK_CHIPX;
         break;

         case HDET_CHIP_Y:
            mask |= HDET_MASK_CHIPY;
         break;

         case HDET_TDET_X:
            mask |= HDET_MASK_TDETX;
         break;

         case HDET_TDET_Y:
            mask |= HDET_MASK_TDETY;
         break;

         case HDET_CHIP_ID:
            mask |= HDET_MASK_CHIPID; 
         break;

         default:
            /* do nothing */ 
         break; 
      } 
   }

   return (mask); 

}


/*H***********************************************************************
 
* DESCRIPTION: The routine dependency_check_hrc is called upon by 
  hrc_process_events to determine if all of the data dependencies required 
  for the requested transformations have been met. The routine checks the 
  data dependencies that need to be verified based on the start of the 
  coordinate transformations. The routine does not validate the values of
  data fields, but merely checks for the existence of the data columns in 
  the input file. A masked value is passed back to the calling routine to 
  be bitwise OR'd to the calling routines error status mask. Missing data 
  dependencies are recorded in the error message stack.

**************************************************************************/

void dependency_check_hrc(
   INPUT_PARMS_P_T inp_p,  /* I - ptr to input parameters data structure */
   STATISTICS_P_T  stat_p, /* I - ptr to statistics data structure       */ 
   dsErrList*      err_p)  /* O - ptr to error message stack             */ 
{
   if (inp_p->stop != HDET_NONE_VAL)
   {
      switch (inp_p->start)
      {
         case HDET_TDET_VAL:
            if ((stat_p->dependencies & HDET_MASK_TDET_REQ) !=
                HDET_MASK_TDET_REQ)
            {
               dsErrAdd(err_p, dsHPEDEPENDENCYERR, Individual, Generic,
                  "TDET requirements");
            }
         break;

         case HDET_CHIP_VAL:
            if (inp_p->hrc_system == HRC_IMAGE_INST)
            {
               if ((stat_p->dependencies & HDET_MASK_CHIP_REQ) !=
                   HDET_MASK_CHIP_REQ)
               {
                  dsErrAdd(err_p, dsHPEDEPENDENCYERR, Individual, Generic,
                     "chipx or chipy columns");
               }
            }
            else
            {
               if ((stat_p->dependencies & HDET_MASK_SCHIP_REQ) !=
                   HDET_MASK_SCHIP_REQ)
               {
                  dsErrAdd(err_p, dsHPEDEPENDENCYERR, Individual, Generic,
                     "chipx/chipy or chipid columns");
               }
            }
         break;

         case HDET_COARSE_VAL:
            if ((stat_p->dependencies & HDET_MASK_CRS_REQ) !=
                HDET_MASK_CRS_REQ)
            {
               dsErrAdd(err_p, dsHPEDEPENDENCYERR, Individual, Generic,
                  "au1-3/av1-3, crsu or crsv columns");
            }
         break;

         default:
            dsErrAdd(err_p, dsHPEDEPENDENCYERR, Individual, Custom,
               "ERROR: The input file is missing data needed to process it in the manner specified.");
         break;
      }
   }
}


/*H***********************************************************************
 
* DESCRIPTION: The routine output_coord_validate is called upon by 
  hrc_process_events to ensure that uncalculated detector coordinate 
  columns are not echoed to the output file. The routine validates the 
  output eventdef list with the coordinate transformation flags and 
  returns a value of TRUE if an invalid output column is listed in the 
  eventdef list.

**************************************************************************/

boolean output_coord_validate(
   short* out_evt_map,  /* I - event column name output position mapping */
   int    out_num,      /* I - number of columns of data in out_evt_map  */
   INPUT_PARMS_P_T inp_p) /* I - data structure containing input params  */  
{
   boolean error = FALSE;
   short rr;
 
   for (rr = 0; rr < out_num; rr++)
   {
      if (((out_evt_map[rr] == HDET_SKY_X) ||
          (out_evt_map[rr] == HDET_SKY_Y)) && (inp_p->stop < HDET_SKY_VAL))
      {
         error = TRUE;
      }
      else
      {
         if (((out_evt_map[rr] == HDET_DET_X) ||
             (out_evt_map[rr] == HDET_DET_Y)) && (inp_p->stop < HDET_DET_VAL))
         {
            error = TRUE;
         }
         else
         {
            if (((out_evt_map[rr] == HDET_TDET_X) ||
                (out_evt_map[rr] == HDET_TDET_Y)) && 
                 (inp_p->stop < HDET_TDET_VAL))
            {
               error = TRUE;
            }
         }
      }
   }
 
   return (error);
}



