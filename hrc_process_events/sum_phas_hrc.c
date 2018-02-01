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
 
* FILE NAME: sum_phas_hrc.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: The routine sum_phas_hrc is called upon by hrc_process_events to 
  calculate the sum of the charge amps for a given event. This is merely
  the total of all six amps (three per each of the two planes).  
 
* NOTES:

  wmclaugh@cfa	Mar 28, 1996  First Version.

* REVISION HISTORY:
  JCC(5/11/00)- compute sum_tot from amps_dd whose data type is in 'double'
*H**************************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 



void sum_phas_hrc(EVENT_REC_T *evt_p) /* I/O structure holding event data */
{
   short count;

   for (count = HDET_1ST_AMP; count < HDET_NUM_AMPS; count++)
   {
      evt_p->amp_tot[HDET_PLANE_X] = evt_p->amp_tot[HDET_PLANE_X] + 
                             evt_p->amps_dd[HDET_PLANE_X][count];  

      evt_p->amp_tot[HDET_PLANE_Y] = evt_p->amp_tot[HDET_PLANE_Y] + 
                             evt_p->amps_dd[HDET_PLANE_Y][count];  
   }

   evt_p->sum_amps = (unsigned short) (evt_p->amp_tot[HDET_PLANE_X] +
                             evt_p->amp_tot[HDET_PLANE_Y]);
}

