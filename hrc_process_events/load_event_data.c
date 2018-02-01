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

/*H**************************************************************************
 
* FILE NAME: load_event_data.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: Load_event_data interates through the evtin_p->num_cols of an event from
  an input qpoe file and copies data into the appropriate fields of an event
  record structure (EVENT_REC_T). The mapping of data fields read in from the
  event to the event record structure are based on a mapping defined in the 
  parameter evtin_p->mapping. This routine is called by hrc_process_events to load 
  the event record structure which will be passed to other routines to calculate
  coordinate values.  
 
* NOTES:

  wmclaugh@cfa	Mar 28, 1996  First Version.

* REVISION HISTORY:
  JCC(5/2/00) - updated to initialize amps_3RD_raw as amps[u:v][HDET_3RD_AMP]
  JCC(5/11/00)- initialize amps_dd as amps_sh ; amps_dd for computation; 
                amps_sh for output;
  JCC(3/26/01)- bugfixed to the Status column which was incorrect when the 
                infile (evt0) had STATUS instead of QUALITY column (eg. 
                outfile from hrc_correct_time).
  JCC(7/2003)-add a new function 'initial_status'
10/2009 - fix rawpos (see Note)
*H**************************************************************************/

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 


void load_event_data(
   EVENT_SETUP_P_T evtin_p,    /* I - input file information               */
   EVENT_REC_P_T  evt_p)       /* I/O structure holding event data         */
{
   int count = 0; 
   boolean vstat = FALSE,  
           etrig = FALSE; 

   while (count < evtin_p->num_cols)
   {
      switch (evtin_p->mapping[count])    
      {
         case HDET_TIME:
            evt_p->time = dmGetScalar_d(evtin_p->desc[count]); 
         break;

         case HDET_MJR_FRAME: 
            evt_p->major_frame = dmGetScalar_l(evtin_p->desc[count]); 
         break;

         case HDET_MNR_FRAME: 
            evt_p->minor_frame = dmGetScalar_l(evtin_p->desc[count]); 
         break;

         case HDET_EVENT: 
            evt_p->event = dmGetScalar_s(evtin_p->desc[count]); 
         break;

         case HDET_CP_X: 
            evt_p->cp[HDET_PLANE_X] = dmGetScalar_s(evtin_p->desc[count]); 
         break;

         case HDET_CP_Y: 
            evt_p->cp[HDET_PLANE_Y] = dmGetScalar_s(evtin_p->desc[count]); 
         break;

         /* JCC(5/11/00) - begin */
         case HDET_AX_1: 
            evt_p->amps_sh[HDET_PLANE_X][HDET_1ST_AMP] = 
               dmGetScalar_s(evtin_p->desc[count]); 
            evt_p->amps_dd[HDET_PLANE_X][HDET_1ST_AMP] = 
               (double)(evt_p->amps_sh[HDET_PLANE_X][HDET_1ST_AMP]) ;
         break;

         case HDET_AX_2: 
            evt_p->amps_sh[HDET_PLANE_X][HDET_2ND_AMP] = 
               dmGetScalar_s(evtin_p->desc[count]); 
            evt_p->amps_dd[HDET_PLANE_X][HDET_2ND_AMP] = 
               (double)(evt_p->amps_sh[HDET_PLANE_X][HDET_2ND_AMP]) ;
         break;
         
         case HDET_AX_3: 
            evt_p->amps_sh[HDET_PLANE_X][HDET_3RD_AMP] = 
               dmGetScalar_s(evtin_p->desc[count]); 
            /* evt_p->amps_3RD_raw[HDET_PLANE_X] = evt_p->amps[HDET_PLANE_X][HDET_3RD_AMP] ;*/
            evt_p->amps_dd[HDET_PLANE_X][HDET_3RD_AMP] = 
               (double)(evt_p->amps_sh[HDET_PLANE_X][HDET_3RD_AMP]) ; 
         break;
         
         case HDET_AY_1: 
            evt_p->amps_sh[HDET_PLANE_Y][HDET_1ST_AMP] = 
               dmGetScalar_s(evtin_p->desc[count]); 
            evt_p->amps_dd[HDET_PLANE_Y][HDET_1ST_AMP] =
               (double)(evt_p->amps_sh[HDET_PLANE_Y][HDET_1ST_AMP]) ; 
         break;
         
         case HDET_AY_2: 
            evt_p->amps_sh[HDET_PLANE_Y][HDET_2ND_AMP] = 
               dmGetScalar_s(evtin_p->desc[count]); 
            evt_p->amps_dd[HDET_PLANE_Y][HDET_2ND_AMP] =
               (double)(evt_p->amps_sh[HDET_PLANE_Y][HDET_2ND_AMP]) ; 
         break;

         case HDET_AY_3: 
            evt_p->amps_sh[HDET_PLANE_Y][HDET_3RD_AMP] = 
               dmGetScalar_s(evtin_p->desc[count]); 
            /* evt_p->amps_3RD_raw[HDET_PLANE_Y] = evt_p->amps[HDET_PLANE_Y][HDET_3RD_AMP] ; */
            evt_p->amps_dd[HDET_PLANE_Y][HDET_3RD_AMP] =
               (double)(evt_p->amps_sh[HDET_PLANE_Y][HDET_3RD_AMP]) ;
         break;
         /* JCC(5/11/00) - end   */

         case HDET_PHA: 
            evt_p->pha = dmGetScalar_s(evtin_p->desc[count]); 
         break;

         case HDET_VETO_STATUS: 
            if (evtin_p->types[count] == dmBIT)
            {
               dmGetScalar_bit(evtin_p->desc[count], &evt_p->veto_status);
            }
            else
            {
               evt_p->veto_status = dmGetScalar_s(evtin_p->desc[count]); 
            } 

            if (!vstat)
            {
               short lld = evt_p->veto_status & 0x13;
 
               /* status bit 6 = vstat bit 5 */ 
               lld += (evt_p->veto_status & 0x20) << 1; 

               /* status bit 3 = vstat bit 2 */ 
               lld += (evt_p->veto_status & 0x04) << 1; 

               evt_p->status = lld << 16; 

               vstat = TRUE; 
            } 
         break;

         case HDET_EVENT_STATUS: 
            evt_p->event_status = dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         /*----------------------------------------------------------*/
         /* JCC(3/26/01) - Add .OR. (|=) to initialize 'status',     */
         /*  so we won't lose the infomation on VETOSTT and E_TRIG   */
         /*  which might be already set earlier in 'evt_p->status'.  */
         /*----------------------------------------------------------*/
         case HDET_STATUS:         /* load the STATUS column from infile */ 
            if (evtin_p->types[count] == dmLONG)
            {  
               /* add  "|"   */
               evt_p->status |= dmGetScalar_l(evtin_p->desc[count]);
            } 
            else
            {
               unsigned char val[4];
               dmGetArray_bit(evtin_p->desc[count], val, 4);

               /* add "|"    */
               evt_p->status |= ( val[3] | (val[2] << 8) |         
                  (val[1]  << 16) | (val[0] << 24)  );
            }

         break; 

         case HDET_X_POS:
            evt_p->xpos = dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_Y_POS:
            evt_p->ypos = dmGetScalar_s(evtin_p->desc[count]); 
         break; 
 
         case HDET_SUMAMPS:
            evt_p->sum_amps = dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_CHIP_ID:
            evt_p->chipid = dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_DUMMY:
            evt_p->dummy = dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_PI:
            evt_p->pi = dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_PHASCALE:
            evt_p->phascale = dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_RAWPHA:
            evt_p->rawpha = dmGetScalar_s(evtin_p->desc[count]); 
         break; 


/* 10/2009-Replace dmGetScalar_s w/ dmGetScalar_l to read evt's 
 *         columns "rawpos(rawx,rawy)"  (ie. HDET_RAW_X, HDET_RAW_Y)
 * Notes :
 * - evt_p->rawpos[0-1] from below are NOT really used anywhere in the code. 
 *   Instead, it always uses the ones that are computed from l1h_coarse_to_raw().
 *   Therefore, this fix should NOT affect any hpe's output.
 */
         case HDET_RAW_X:
evt_p->rawpos[HDET_PLANE_X]=dmGetScalar_l(evtin_p->desc[count]);/*10/2009-change _s to _l*/
         break; 

         case HDET_RAW_Y:
evt_p->rawpos[HDET_PLANE_Y]=dmGetScalar_l(evtin_p->desc[count]);/*10/2009-change _s to _l*/
         break; 
/*end:*/


         case HDET_DET_X:
            evt_p->detpos[HDET_PLANE_X] = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_DET_Y:
            evt_p->detpos[HDET_PLANE_Y] = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_TDET_X:
            evt_p->tdetpos[HDET_PLANE_X] = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_TDET_Y:
            evt_p->tdetpos[HDET_PLANE_Y] = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_SKY_X:
            evt_p->skypos[HDET_PLANE_X] =
               dmGetScalar_s(evtin_p->desc[count]);
         break;
 
         case HDET_SKY_Y:
            evt_p->skypos[HDET_PLANE_Y] =
               dmGetScalar_s(evtin_p->desc[count]);
         break;

         case HDET_SKY:
            dmGetVector_d(evtin_p->desc[count], evt_p->skypos, 2);
         break; 

         case HDET_FPZ:
            evt_p->skypos[HDET_PLANE_Z] = 
               dmGetScalar_l(evtin_p->desc[count]); 
         break; 

         case HDET_CHIP_X:
            evt_p->chippos[HDET_PLANE_X] = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_CHIP_Y:
            evt_p->chippos[HDET_PLANE_Y] = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 
   
         case HDET_TICK:
            evt_p->tick = 
               dmGetScalar_l(evtin_p->desc[count]); 
         break; 
   
         case HDET_SCIFR:
            evt_p->scifr = 
               dmGetScalar_l(evtin_p->desc[count]); 
         break; 
   
         case HDET_EVTCTR:
            evt_p->evtctr = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 
   
         case HDET_AMP_SF:
            evt_p->amp_sf = 
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_VETO_STT:
            if (evtin_p->types[count] == dmBIT)
            {
               dmGetScalar_bit(evtin_p->desc[count], &evt_p->veto_stt);
            }
            else
            {
               evt_p->veto_stt = 
                  dmGetScalar_s(evtin_p->desc[count]); 
            } 
            if (!vstat)
            {
               short lld = evt_p->veto_stt & 0xFE;

               /* negate bit 0 */ 
               lld += !(evt_p->veto_stt & 0x01); 

               evt_p->status = lld << 16;  

               vstat = TRUE;
            } 
         break; 

         case HDET_E_TRIG:
            if (evtin_p->types[count] == dmBIT)
            {
               dmGetScalar_bit(evtin_p->desc[count], &evt_p->e_trig);
            }
            else 
            {
            
               evt_p->e_trig = dmGetScalar_s(evtin_p->desc[count]); 
            } 
            etrig = TRUE; 
         break; 
   
         case HDET_SUBMJF:
            evt_p->submjf =
               dmGetScalar_s(evtin_p->desc[count]); 
         break; 

         case HDET_DET_ID:
            if (evtin_p->types[count] == dmBIT)
            {
               dmGetScalar_bit(evtin_p->desc[count], &evt_p->det_id);
            } 
            else
            {
               evt_p->det_id = dmGetScalar_s(evtin_p->desc[count]); 
            } 
         break;

         case HDET_MRF:
            evt_p->mrf =
               dmGetScalar_s(evtin_p->desc[count]); 
         break;

         case HDET_STOPMNF:
            evt_p->stopmnf =
               dmGetScalar_s(evtin_p->desc[count]); 
         break;

   
         default:
            /* do nothing- we don't care about extraneous data fields */
            /*    in input data since we only need to output desired  */
            /*    fields.                                             */  
         break; 

      } /* end switch */ 

      count++; 
   } /* end while */ 

   if (etrig)
   {
      /* negate bits 1,0 from e_trig */ 
      evt_p->status |= (!(evt_p->e_trig & 0x02)) << 25;     
      evt_p->status |= (!(evt_p->e_trig & 0x01)) << 24;     
   }
   else
   {
      /* negate event_status bit 7 and veto_stat bit 6 */ 
      evt_p->status |= (evt_p->event_status & 0x80) << 18; 
      evt_p->status |= (evt_p->veto_status & 0x40) << 18; 
   } 
} /* load_event_data */


/*******************************************************************
 *JCC(7/2003)
 * -initalize status bits 6 to 15 and 26 
 *******************************************************************/
void initial_status(
   INPUT_PARMS_P_T  inp_p,     /* I: do_ratio */
   EVENT_REC_P_T  evt_p)       /* I/O structure holding event data */
{
   evt_p->status  &= ~(HDET_SEQUENCE_STS | HDET_PI_VALUE_STS |
                       HDET_ZERO_SUM_STS | HDET_U_CNTR_STS |
                       HDET_V_CNTR_STS | HDET_FIN_POS_STS | 
                       HDET_HOT_SPOT_STS ) ;

   if (inp_p->do_ratio)        /* true */
   {
      evt_p->status  &= ~( HDET_GRID_RTO_STS | HDET_ZERO_PSUM_STS |
                           HDET_PHA_RTO_STS ) ;
   }
} /* initial_status */
