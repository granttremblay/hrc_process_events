/*                                                                
**  Copyright (C) 1996-2009  Smithsonian Astrophysical Observatory 
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
 
* FILE NAME: write_hrc_events.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: The routine write_hrc_events is called by hrc_process_events
  to write a specified event to an output qpoe file. The attributes which
  are to be output are specified in the hrc_process_events.par paramater
  file. A mapping array is passed into this routine to facilitate the
  output of the appropriate data fields. The routine does not return any
  error status but adds detected errors onto the error list which is passed
  in.
 
* NOTES:
 
  wmclaugh@cfa  Mar 28, 1996  First Version.
 
* REVISION HISTORY:
  JCC(5/1/00)-the raw A3 is now stored in amps_3RD_raw.
  JCC(5/11/00)- write the original raw amplitudes (amps_sh)
  JCC(6/18/00)- add comment for HDET_STATUS .
  JCC(8/2/00)-add FLOAT for det,sky coords.
10/2009-dph/fap new gain files affect PI (see 'Notes on outCol PI')
*H***********************************************************************/

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 


void write_hrc_events(
   EVENT_SETUP_P_T evtout_p, /* I - output event file info              */
   EVENT_REC_P_T   evt_p,    /* I - structure holding event data        */
   dsErrList*      err_p)    /* O - error list pointer                  */
{
   int count = 0;

   while (count < evtout_p->num_cols)
   {
      switch (evtout_p->mapping[count])
      {
         case HDET_TIME:
            dmSetScalar_d(evtout_p->desc[count], evt_p->time); 
         break;

         case HDET_MJR_FRAME:
            dmSetScalar_l(evtout_p->desc[count], evt_p->major_frame); 
         break;

         case HDET_MNR_FRAME:
            dmSetScalar_l(evtout_p->desc[count], evt_p->minor_frame); 
         break;

         case HDET_EVENT:
            dmSetScalar_s(evtout_p->desc[count], evt_p->event); 
         break;

         case HDET_CP_X:
            dmSetScalar_s(evtout_p->desc[count], evt_p->cp[HDET_PLANE_X]); 
         break;

         case HDET_CP_Y:
            dmSetScalar_s(evtout_p->desc[count], evt_p->cp[HDET_PLANE_Y]); 
         break;

         /* JCC (5/11/00) - begin */
         case HDET_AX_1:
            dmSetScalar_s(evtout_p->desc[count], 
               evt_p->amps_sh[HDET_PLANE_X][HDET_1ST_AMP]); 
         break;

         case HDET_AX_2:
            dmSetScalar_s(evtout_p->desc[count], 
               evt_p->amps_sh[HDET_PLANE_X][HDET_2ND_AMP]); 
         break;

         case HDET_AX_3:
            dmSetScalar_s(evtout_p->desc[count],
               evt_p->amps_sh[HDET_PLANE_X][HDET_3RD_AMP]); 
 /* dmSetScalar_s(evtout_p->desc[count], evt_p->amps_3RD_raw[HDET_PLANE_X]); */
         break;

         case HDET_AY_1:
            dmSetScalar_s(evtout_p->desc[count], 
               evt_p->amps_sh[HDET_PLANE_Y][HDET_1ST_AMP]); 
         break;

         case HDET_AY_2:
            dmSetScalar_s(evtout_p->desc[count], 
               evt_p->amps_sh[HDET_PLANE_Y][HDET_2ND_AMP]); 
         break;

         case HDET_AY_3:
            dmSetScalar_s(evtout_p->desc[count], 
               evt_p->amps_sh[HDET_PLANE_Y][HDET_3RD_AMP]); 
 /* dmSetScalar_s( evtout_p->desc[count], evt_p->amps_3RD_raw[HDET_PLANE_Y]); */
         break; 
         /* JCC (5/11/00) - end */

         case HDET_PHA:
            dmSetScalar_s(evtout_p->desc[count], evt_p->pha); 
         break;

         case HDET_VETO_STATUS:
            if (evtout_p->types[count] == dmBIT)
            {
               dmSetArray_bit(evtout_p->desc[count], &evt_p->veto_status, 1);
            }
            else
            {
               dmSetScalar_s(evtout_p->desc[count], evt_p->veto_status); 
            }
         break;

         case HDET_EVENT_STATUS:
            dmSetScalar_s(evtout_p->desc[count], evt_p->event_status); 
         break;

         case HDET_STATUS:        /*output column 'status' */ 
            if (evtout_p->types[count] == dmBIT) 
            {
               /*----------------------------------------------------------
                *JCC(6/16/00)-data type of 'status' is 'bit' (ie. x:status)
                * 
                * val[0] = (evt_p->status & 0xff000000) >> 24;
                *  This takes the high end value(0xff000000) from status
                *  and shift it to the right for 3 bytes (24/8=3).
                *---------------------------------------------------------*/
               unsigned char val[4];

               /* to avoid byte swapping on little endian platforms */ 
               val[0] = (evt_p->status & 0xff000000) >> 24;  
               val[1] = (evt_p->status & 0x00ff0000) >> 16;  
               val[2] = (evt_p->status & 0x0000ff00) >> 8;  
               val[3] = evt_p->status & 0x000000ff;  

               dmSetArray_bit(evtout_p->desc[count], val, 4); 
            } 
            else
            {
               /*----------------------------------------------------------
                *JCC(6/16/00)-data type of 'status' is 'long' (ie. l:status)
                *---------------------------------------------------------*/
               dmSetScalar_l(evtout_p->desc[count], evt_p->status); 
            } 
         break; 

         case HDET_X_POS:
            dmSetScalar_s(evtout_p->desc[count], evt_p->xpos); 
         break;

         case HDET_Y_POS:
            dmSetScalar_s(evtout_p->desc[count], evt_p->ypos); 
         break;

         case HDET_SUMAMPS:
            dmSetScalar_s(evtout_p->desc[count], evt_p->sum_amps); 
         break;

         case HDET_CHIP_ID:
            dmSetScalar_s(evtout_p->desc[count], evt_p->chipid); 
         break;

         case HDET_DUMMY:
            dmSetScalar_s(evtout_p->desc[count], evt_p->dummy); 
         break;

         case HDET_PI: /*11/2009: outCol maxPI=255||1023; SHORT dtype*/
            dmSetScalar_s(evtout_p->desc[count], (short)evt_p->pi);
         break;

         case HDET_PHASCALE:
            dmSetScalar_s(evtout_p->desc[count], evt_p->phascale); 
         break;

         case HDET_RAWPHA:
            dmSetScalar_s(evtout_p->desc[count], evt_p->rawpha); 
         break;

         case HDET_TICK:
            dmSetScalar_l(evtout_p->desc[count], evt_p->tick); 
         break;

         case HDET_SCIFR:
            dmSetScalar_l(evtout_p->desc[count], evt_p->scifr); 
         break;

         case HDET_EVTCTR:
            dmSetScalar_s(evtout_p->desc[count], evt_p->evtctr); 
         break;

         case HDET_RAW_X:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count],
                                evt_p->rawpos[HDET_PLANE_X]);
               break;
 
               case dmLONG:
               {
                  long cast_val;
 
                  pix_double_to_long(evt_p->rawpos[HDET_PLANE_X],
                                     FALSE, &cast_val);
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break;
 
               case dmSHORT: /* fallthrough intented */
               default:
               {
                  short cast_val;
 
                  pix_double_to_short(evt_p->rawpos[HDET_PLANE_X],
                                      FALSE, &cast_val);
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break;
            }
         break;

         case HDET_RAW_Y:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count],
                                evt_p->rawpos[HDET_PLANE_Y]);
               break;
 
               case dmLONG:
               {
                  long cast_val;
 
                  pix_double_to_long(evt_p->rawpos[HDET_PLANE_Y],
                                     FALSE, &cast_val);
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break;
 
               case dmSHORT: /* fallthrough intented */
               default:
               {
                  short cast_val;
 
                  pix_double_to_short(evt_p->rawpos[HDET_PLANE_Y],
                                      FALSE, &cast_val);
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break;
            }
         break;

         case HDET_RAW:
            switch (evtout_p->types[count])
            {
               case dmSHORT:
               {
                  short casted_val[2];
 
                  pix_double_to_short(evt_p->rawpos[HDET_PLANE_X],
                                      FALSE, &casted_val[0]);
                  pix_double_to_short(evt_p->rawpos[HDET_PLANE_Y],
                                      FALSE, &casted_val[1]);
                  dmSetVector_s(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmLONG:
               {
                  long casted_val[2];
 
                  pix_double_to_long(evt_p->rawpos[HDET_PLANE_X],
                                     FALSE, &casted_val[0]);
                  pix_double_to_long(evt_p->rawpos[HDET_PLANE_Y],
                                     FALSE, &casted_val[1]);
                  dmSetVector_l(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmDOUBLE:
                  dmSetVector_d(evtout_p->desc[count], evt_p->rawpos, 2);
               break;
 
               default:
               break;
            } /* end switch */
         break; 

         case HDET_DET_X:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count],
                                evt_p->detpos[HDET_PLANE_X]);
               break;

               case dmFLOAT:
               {
                  float cast_val;

                  cast_val = (float)(evt_p->detpos[HDET_PLANE_X]) ;
                  dmSetScalar_f(evtout_p->desc[count], cast_val);
               }
               break;

               case dmLONG:
               {
                  long cast_val; 

                  pix_double_to_long(evt_p->detpos[HDET_PLANE_X], 
                                     FALSE, &cast_val);
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break; 
 
               case dmSHORT: /* fallthrough intented */
               default:
               {
                  short cast_val; 

                  pix_double_to_short(evt_p->detpos[HDET_PLANE_X], 
                                      FALSE, &cast_val);
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break;
            }
         break;

         case HDET_DET_Y:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count],
                                evt_p->detpos[HDET_PLANE_Y]);
               break;

               case dmFLOAT:
               {
                  float cast_val;

                  cast_val = (float)(evt_p->detpos[HDET_PLANE_Y]);
                  dmSetScalar_f(evtout_p->desc[count], cast_val);
               }
               break;

               case dmLONG:
               {
                  long cast_val; 

                  pix_double_to_long(evt_p->detpos[HDET_PLANE_Y], 
                                      FALSE, &cast_val);
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break; 
 
               case dmSHORT: /* fallthrough intented */
               default:
               {
                  short cast_val; 

                  pix_double_to_short(evt_p->detpos[HDET_PLANE_Y], 
                                      FALSE, &cast_val);
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break;
            }
         break;

         case HDET_DET:
         {
            switch (evtout_p->types[count])
            {
               case dmSHORT:
               {
                  short casted_val[2];
 
                  pix_double_to_short(evt_p->detpos[HDET_PLANE_X],
                                      FALSE, &casted_val[0]);
                  pix_double_to_short(evt_p->detpos[HDET_PLANE_Y],
                                      FALSE, &casted_val[1]);
                  dmSetVector_s(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmLONG:
               {
                  long casted_val[2];
 
                  pix_double_to_long(evt_p->detpos[HDET_PLANE_X],
                                     FALSE, &casted_val[0]);
                  pix_double_to_long(evt_p->detpos[HDET_PLANE_Y],
                                     FALSE, &casted_val[1]);
                  dmSetVector_l(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmDOUBLE:
                  dmSetVector_d(evtout_p->desc[count], evt_p->detpos, 2);
               break;
 
               case dmFLOAT:
               {
                  float casted_val[2];

                  casted_val[0] = (float)(evt_p->detpos[HDET_PLANE_X]);
                  casted_val[1] = (float)(evt_p->detpos[HDET_PLANE_Y]);
                  dmSetVector_f(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               default:
               break;
            } /* end switch */
         }
         break;

         case HDET_TDET_X:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count],
                                evt_p->tdetpos[HDET_PLANE_X]);
               break;

               case dmLONG:
               {
                  long cast_val;

                  pix_double_to_long(evt_p->tdetpos[HDET_PLANE_X], FALSE,
                                     &cast_val);  
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break; 
 
               case dmSHORT: /* fallthrough intented */
               default:
               {
                  short cast_val;

                  pix_double_to_short(evt_p->tdetpos[HDET_PLANE_X], FALSE,
                                      &cast_val);  
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break;
            }
         break;
 
         case HDET_TDET_Y:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count],
                                evt_p->tdetpos[HDET_PLANE_Y]);
               break;

               case dmLONG:
               {
                  long cast_val; 

                  pix_double_to_long(evt_p->tdetpos[HDET_PLANE_Y], FALSE,
                                     &cast_val);  
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break; 
 
               case dmSHORT: /* fallthrough intented */
               default:
               {
                  short cast_val; 

                  pix_double_to_short(evt_p->tdetpos[HDET_PLANE_Y], FALSE,
                                      &cast_val);  
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break;
            }
         break;

         case HDET_TDET:
         {
            switch (evtout_p->types[count])
            {
               case dmSHORT:
               {
                  short casted_val[2];
 
                  pix_double_to_short(evt_p->tdetpos[HDET_PLANE_X],
                                      FALSE, &casted_val[0]);
                  pix_double_to_short(evt_p->tdetpos[HDET_PLANE_Y],
                                      FALSE, &casted_val[1]);
                  dmSetVector_s(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmLONG:
               {
                  long casted_val[2];

                  pix_double_to_long(evt_p->tdetpos[HDET_PLANE_X],
                                     FALSE, &casted_val[0]);
                  pix_double_to_long(evt_p->tdetpos[HDET_PLANE_Y],
                                     FALSE, &casted_val[1]);
                  dmSetVector_l(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmDOUBLE:
                  dmSetVector_d(evtout_p->desc[count], evt_p->tdetpos, 2);
               break;
 
               default:
               break;
            } /* end switch */
         }
         break;

         case HDET_SKY_X:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count], 
                                evt_p->skypos[HDET_PLANE_X]);
               break;

               case dmFLOAT:
               {
                  float cast_val;

                  cast_val = (float)(evt_p->skypos[HDET_PLANE_X]);
                  dmSetScalar_f(evtout_p->desc[count], cast_val);
               }
               break;

               case dmLONG:  
               {
                  long cast_val;

                  pix_double_to_long(evt_p->skypos[HDET_PLANE_X],
                                     FALSE, &cast_val);
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break; 

               case dmSHORT: /* fallthrough intented */ 
               default: 
               {
                  short cast_val;

                  pix_double_to_short(evt_p->skypos[HDET_PLANE_X], 
                                      FALSE, &cast_val);
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break; 
            }
         break;
 
         case HDET_SKY_Y:
            switch(evtout_p->types[count])
            {
               case dmDOUBLE:
                  dmSetScalar_d(evtout_p->desc[count], 
                                evt_p->skypos[HDET_PLANE_Y]);
               break;

               case dmFLOAT:
               {
                  float cast_val;

                  cast_val = (float)(evt_p->skypos[HDET_PLANE_Y]);
                  dmSetScalar_f(evtout_p->desc[count], cast_val);
               }
               break;

               case dmLONG:
               {
                  long cast_val;

                  pix_double_to_long(evt_p->skypos[HDET_PLANE_Y], 
                                     FALSE, &cast_val);
#ifdef AXLEN_WORKAROUND
                  cast_val -= 28671; 
#endif 
                  dmSetScalar_l(evtout_p->desc[count], cast_val);
               }
               break; 

               case dmSHORT: /* fallthrough intented */
               default:
               {
                  short cast_val;

                  pix_double_to_short(evt_p->skypos[HDET_PLANE_Y], 
                                         FALSE, &cast_val);
                  dmSetScalar_s(evtout_p->desc[count], cast_val);
               }
               break;
            }
         break;

         case HDET_SKY:
         {
            switch (evtout_p->types[count])
            {
               case dmSHORT:
               {
                  short casted_val[2];
 
                  casted_val[0] = (float)(evt_p->skypos[HDET_PLANE_X]);
                  casted_val[1] = (float)(evt_p->skypos[HDET_PLANE_Y]);
                  dmSetVector_s(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmLONG:
               {
                  long casted_val[2];
 
                  pix_double_to_long(evt_p->skypos[HDET_PLANE_X],
                                     FALSE, &casted_val[0]);
                  pix_double_to_long(evt_p->skypos[HDET_PLANE_Y],
                                     FALSE, &casted_val[1]);
                  dmSetVector_l(evtout_p->desc[count], casted_val, 2);
               }
               break;
 
               case dmDOUBLE:
                  dmSetVector_d(evtout_p->desc[count], evt_p->skypos, 2);
               break;
 
               case dmFLOAT:
               {
                  float casted_val[2];

                  casted_val[0] = (float)(evt_p->skypos[HDET_PLANE_X]);
                  casted_val[1] = (float)(evt_p->skypos[HDET_PLANE_Y]);
                  dmSetVector_f(evtout_p->desc[count], casted_val, 2); 
               }
               break;
 
               default:
               break;
            } /* end switch */
         }
         break;

         case HDET_FPZ:
         {
            long cast_val;

            pix_double_to_long(evt_p->skypos[HDET_PLANE_Z], 
                               FALSE, &cast_val);
            dmSetScalar_l(evtout_p->desc[count], cast_val); 
         } 
         break; 

         case HDET_CHIP_X:
         {
            short cast_val; 

            pix_double_to_short(evt_p->chippos[HDET_PLANE_X], FALSE,
                                &cast_val); 
            dmSetScalar_s(evtout_p->desc[count], cast_val);
         }
         break;

         case HDET_CHIP_Y:
         {
            short cast_val; 

            pix_double_to_short(evt_p->chippos[HDET_PLANE_Y], FALSE,
                                &cast_val); 
            dmSetScalar_s(evtout_p->desc[count], cast_val);
         }
         break;

         case HDET_CHIP:
         {
            short casted_val[2];
 
            pix_double_to_short(evt_p->chippos[HDET_PLANE_X],
                                FALSE, &casted_val[0]);
            pix_double_to_short(evt_p->chippos[HDET_PLANE_Y],
                                FALSE, &casted_val[1]);
 
            dmSetVector_s(evtout_p->desc[count], casted_val, 2);
         } 
         break;


         case HDET_AMP_SF:
            dmSetScalar_s(evtout_p->desc[count], evt_p->amp_sf);
         break;

         case HDET_E_TRIG:
            if (evtout_p->types[count] == dmBIT)
            {
               dmSetArray_bit(evtout_p->desc[count], &evt_p->e_trig, 1); 
            } 
            else
            {
               dmSetScalar_s(evtout_p->desc[count], evt_p->e_trig);
            } 

         break; 

         case HDET_VETO_STT:
            if (evtout_p->types[count] == dmBIT)
            {
               dmSetArray_bit(evtout_p->desc[count], &evt_p->veto_stt, 1);
            }
            else
            {
               dmSetScalar_s(evtout_p->desc[count], evt_p->veto_stt);
            } 
         break;

         case HDET_SUBMJF:
            dmSetScalar_s(evtout_p->desc[count], evt_p->submjf);
         break; 

         case HDET_DET_ID:
            if (evtout_p->types[count] == dmBIT)
            {
               dmSetArray_bit(evtout_p->desc[count], &evt_p->det_id, 1);
            }
            else
            {
               dmSetScalar_s(evtout_p->desc[count], evt_p->det_id);
            } 
         break; 

         case HDET_STOPMNF:
            dmSetScalar_s(evtout_p->desc[count], evt_p->stopmnf);
         break; 

         case HDET_MRF:
            dmSetScalar_s(evtout_p->desc[count], evt_p->mrf);
         break; 

         default:
            /* unknown data field for output */
            dsErrAdd(err_p, dsHPEWRITEEVTERR, Accumulation, Generic,
               evtout_p->file);
         break; 
      } /* end switch */

      count++;
   } /* end while */

   /* write out event */ 
   dmTablePutRow(evtout_p->extension, NULL);
}

