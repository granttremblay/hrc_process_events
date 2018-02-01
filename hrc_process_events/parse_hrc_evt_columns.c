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

/*H**************************************************************************
 
* FILE NAME: parse_hrc_evt_columns.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: Parse_hrc_evt_columns is called on by hrc_process_events to map 
  given data columns to event record fields. The routine is used to set up both
  the input mapping and the output mapping. Mappings are done based on string
  macros defined in hrc_process_events.h. If greater flexibility is desired in 
  the mappings, new defines can be added to the .h file and extra conditionals
  can be added to this routine. For instance, the read coordinate for x can be
  either "readx" or "read_x".
 
  The routine returns a boolean value of TRUE if an error is encountered (ie.
  it does not recognize a particular field) and a FALSE otherwise. 
 
* NOTES:

  wmclaugh@cfa	Mar 28, 1996  First Version.

* REVISION HISTORY:
 
        Ref. No.        Date
        --------        ----
        1.1             28 Mar 1996
*H**************************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 


bool parse_hrc_evt_columns(
   char** out_att_names,  /* I - event file column names */  
   int    columns,        /* I - number of columns in event */ 
   short* out_evt_map)    /* O - event column name output position mapping */ 
{
   bool error = FALSE; 
   short count = 0;
   short pos; 

   while (count < columns)  
   {
      char  name_uc[DS_SZ_COLNAME];
 
      memset(name_uc, 0, DS_SZ_COLNAME*sizeof(char));

      /* convert name to uppercase before comparing */
      if ((pos = strlen(out_att_names[count])) > DS_SZ_PATHNAME)
      {
         pos = DS_SZ_PATHNAME - 1;
      }
      while (pos-- > 0)
      {
         name_uc[pos] = (char) toupper(out_att_names[count][pos]);
      }

      if ((strcmp(name_uc, HDET_TIME_COL)) == 0) 
      {
         out_evt_map[count] = HDET_TIME; 
      }
      else if ((strcmp(name_uc, HDET_MJR_FRAME_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_MJR_FRAME; 
      } 
      else if ((strcmp(name_uc, HDET_MJR_FRAME_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_MJR_FRAME; 
      } 
      else if ((strcmp(name_uc, HDET_MNR_FRAME_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_MNR_FRAME; 
      } 
      else if ((strcmp(name_uc, HDET_MNR_FRAME_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_MNR_FRAME; 
      } 
      else if ((strcmp(name_uc, HDET_MNR_FRAME_COL_ALT2)) == 0) 
      { 
         out_evt_map[count] = HDET_MNR_FRAME; 
      } 
      else if ((strcmp(name_uc, HDET_MNR_FRAME_COL_ALT3)) == 0) 
      { 
         out_evt_map[count] = HDET_MNR_FRAME; 
      } 
      else if ((strcmp(name_uc, HDET_EVENT_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_EVENT; 
      } 
      else if ((strcmp(name_uc, HDET_CP_X_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_X; 
      } 
      else if ((strcmp(name_uc, HDET_CP_X_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_X; 
      } 
      else if ((strcmp(name_uc, HDET_CP_X_COL_ALT2)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_X; 
      } 
      else if ((strcmp(name_uc, HDET_CP_X_COL_ALT3)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_X; 
      } 
      else if ((strcmp(name_uc, HDET_CP_Y_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_Y; 
      } 
      else if ((strcmp(name_uc, HDET_CP_Y_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_Y; 
      } 
      else if ((strcmp(name_uc, HDET_CP_Y_COL_ALT2)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_Y; 
      } 
      else if ((strcmp(name_uc, HDET_CP_Y_COL_ALT3)) == 0) 
      { 
         out_evt_map[count] = HDET_CP_Y; 
      } 
      else if ((strcmp(name_uc, HDET_AX_1_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_AX_1; 
      } 
      else if ((strcmp(name_uc, HDET_AX_1_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_AX_1; 
      } 
      else if ((strcmp(name_uc, HDET_AX_2_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_AX_2; 
      } 
      else if ((strcmp(name_uc, HDET_AX_2_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_AX_2; 
      } 
      else if ((strcmp(name_uc, HDET_AX_3_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_AX_3; 
      } 
      else if ((strcmp(name_uc, HDET_AX_3_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_AX_3; 
      } 
      else if ((strcmp(name_uc, HDET_AY_1_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_AY_1; 
      } 
      else if ((strcmp(name_uc, HDET_AY_1_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_AY_1; 
      } 
      else if ((strcmp(name_uc, HDET_AY_2_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_AY_2; 
      } 
      else if ((strcmp(name_uc, HDET_AY_2_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_AY_2; 
      } 
      else if ((strcmp(name_uc, HDET_AY_3_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_AY_3; 
      } 
      else if ((strcmp(name_uc, HDET_AY_3_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_AY_3; 
      } 
      else if ((strcmp(name_uc, HDET_PHA_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_PHA; 
      } 
      else if ((strcmp(name_uc, HDET_V_STS_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_VETO_STATUS; 
      } 
      else if ((strcmp(name_uc, HDET_V_STS_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_VETO_STATUS; 
      } 
      else if ((strcmp(name_uc, HDET_E_STS_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_EVENT_STATUS; 
      } 
      else if ((strcmp(name_uc, HDET_STATUS_COL)) == 0)
      {
         out_evt_map[count] = HDET_STATUS;
      }
      else if ((strcmp(name_uc, HDET_STATUS_COL_ALT1)) == 0)
      {
         out_evt_map[count] = HDET_STATUS;
      }
      else if ((strcmp(name_uc, HDET_X_POS_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_X_POS; 
      } 
      else if ((strcmp(name_uc, HDET_Y_POS_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_Y_POS; 
      } 
      else if ((strcmp(name_uc, HDET_SUMAMPS_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_SUMAMPS; 
      } 
      else if ((strcmp(name_uc, HDET_CHIP_ID_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_CHIP_ID; 
      } 
      else if ((strcmp(name_uc, HDET_DUMMY_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_DUMMY; 
      } 
      else if ((strcmp(name_uc, HDET_PI_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_PI; 
      } 
      else if ((strcmp(name_uc, HDET_PHASCALE_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_PHASCALE; 
      } 
      else if ((strcmp(name_uc, HDET_RAWPHA_COL)) == 0) 
      { 
         /* not sure if raw_pha and pha should be the same entry */
         /* this seems to be the case from the sample run with the original */
         /* code, so for now, comment out the raw_pha line and replace with */
         /* the pha line... may need to change later */ 
         /*
         out_evt_map[count] = HDET_RAWPHA; 
         */ 

         out_evt_map[count] = HDET_PHA; 
      } 
      else if ((strcmp(name_uc, HDET_RAW_X_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_RAW_X; 
      } 
      else if ((strcmp(name_uc, HDET_RAW_X_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_RAW_X; 
      } 
      else if ((strcmp(name_uc, HDET_RAW_Y_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_RAW_Y; 
      } 
      else if ((strcmp(name_uc, HDET_RAW_Y_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_RAW_Y; 
      } 
      else if ((strcmp(name_uc, HDET_RAW_COL)) == 0)
      {
         out_evt_map[count] = HDET_RAW;
      }
      else if ((strcmp(name_uc, HDET_DET_X_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_DET_X; 
      } 
      else if ((strcmp(name_uc, HDET_DET_X_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_DET_X; 
      } 
      else if ((strcmp(name_uc, HDET_DET_Y_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_DET_Y; 
      } 
      else if ((strcmp(name_uc, HDET_DET_Y_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_DET_Y; 
      } 
      else if ((strcmp(name_uc, HDET_DET_COL)) == 0)
      {
         out_evt_map[count] = HDET_DET;
      }
      else if ((strcmp(name_uc, HDET_CHIP_X_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_CHIP_X; 
      } 
      else if ((strcmp(name_uc, HDET_CHIP_X_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_CHIP_X; 
      } 
      else if ((strcmp(name_uc, HDET_CHIP_Y_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_CHIP_Y; 
      } 
      else if ((strcmp(name_uc, HDET_CHIP_Y_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_CHIP_Y; 
      } 
      else if ((strcmp(name_uc, HDET_CHIP_COL)) == 0)
      {
         out_evt_map[count] = HDET_CHIP;
      }
      else if ((strcmp(name_uc, HDET_TDET_X_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_TDET_X; 
      } 
      else if ((strcmp(name_uc, HDET_TDET_X_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_TDET_X; 
      } 
      else if ((strcmp(name_uc, HDET_TDET_Y_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_TDET_Y; 
      } 
      else if ((strcmp(name_uc, HDET_TDET_Y_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_TDET_Y; 
      } 
      else if ((strcmp(name_uc, HDET_TDET_COL)) == 0)
      {
         out_evt_map[count] = HDET_TDET;
      }
      else if ((strcmp(name_uc, HDET_SKY_X_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_SKY_X; 
      } 
      else if ((strcmp(name_uc, HDET_SKY_Y_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_SKY_Y; 
      } 
      else if ((strcmp(name_uc, HDET_SKY_COL)) == 0)
      {
         out_evt_map[count] = HDET_SKY;
      }
      else if ((strcmp(name_uc, HDET_FPZ_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_FPZ; 
      } 
      else if ((strcmp(name_uc, HDET_TICK_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_TICK; 
      } 
      else if ((strcmp(name_uc, HDET_TICK_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_TICK; 
      } 
      else if ((strcmp(name_uc, HDET_SCIFR_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_SCIFR; 
      } 
      else if ((strcmp(name_uc, HDET_AMP_SF_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_AMP_SF; 
      } 
      else if ((strcmp(name_uc, HDET_EVTCTR_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_EVTCTR; 
      } 
      else if ((strcmp(name_uc, HDET_E_TRIG_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_E_TRIG; 
      } 
      else if ((strcmp(name_uc, HDET_E_TRIG_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_E_TRIG; 
      } 
      else if ((strcmp(name_uc, HDET_VETO_STT_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_VETO_STT; 
      } 
      else if ((strcmp(name_uc, HDET_VETO_STT_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_VETO_STT; 
      } 
      else if ((strcmp(name_uc, HDET_SUBMJF_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_SUBMJF; 
      } 
      else if ((strcmp(name_uc, HDET_SUBMJF_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_SUBMJF; 
      } 
      else if ((strcmp(name_uc, HDET_DET_ID_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_DET_ID; 
      } 
      else if ((strcmp(name_uc, HDET_MRF_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_MRF; 
      } 
      else if ((strcmp(name_uc, HDET_STOPMNF_COL)) == 0) 
      { 
         out_evt_map[count] = HDET_STOPMNF; 
      } 
      else if ((strcmp(name_uc, HDET_STOPMNF_COL_ALT1)) == 0) 
      { 
         out_evt_map[count] = HDET_STOPMNF; 
      } 
      else
      {
         /* unknown column */ 
         error = TRUE; 
         out_evt_map[count] = HDET_UNKNOWN_FIELD; 
      }
      count++; 
   }

   return (error); 
}

