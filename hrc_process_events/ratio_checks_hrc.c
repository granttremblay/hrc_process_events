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
 
* FILE NAME: ratio_checks_hrc.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: The routine ratio_checks_hrc is called upon by hrc_process_events
  to perform validity checks on the event data. The sum of amps to pha ratio as 
  well as grid charge ratios are checked. If either is invalid, the appropriate
  ratio count is incremented to kep track of the 'bad' event. The checks are 
  only performed if the parameters for each are properly set in the  
  hrc_process_events.par parameter file. 
 
* NOTES:

  The routine makes the assumption that sum_phas_hrc has been called prior to 
  invoking this routine as both ratio checks depend on values computed in
  the sum_phas_hrc routine (namely, sum_phas, amp_tot[HDET_PLANE_X], and
  amp_tot[HDET_PLANE_Y]). 

  wmclaugh@cfa	Mar 28, 1996  First Version.

* REVISION HISTORY:
 
        Ref. No.        Date
        --------        ----
        1.1             28 Mar 1996
*H**************************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 


void ratio_checks_hrc(
   EVENT_REC_T     *evt_p,  /* I   structure holding event data              */
   INPUT_PARMS_P_T  inp_p,  /* I   ptr to struct containing parameter data   */ 
   STATISTICS_T    *stat_p) /* O   structure holding statistical data counts */ 
{
   /* grid charge ratio validity check */ 
   if ((inp_p->grid_ratio != 0) && (evt_p->amp_tot[HDET_PLANE_Y] != 0))
   {
      double ratio = evt_p->amp_tot[HDET_PLANE_X]/evt_p->amp_tot[HDET_PLANE_Y]; 

      if ((ratio > (1 + inp_p->grid_ratio)) ||
          (ratio < (1 - inp_p->grid_ratio))) 
      {
         stat_p->bad_grid_ratio++; 
         evt_p->status |= HDET_GRID_RTO_STS; 
      } 
   }

   /* check sum_amps to pha ratio */
   if (inp_p->pha_ratio != 0)
   {
      int scale; 
      int spha; 
      
      if (inp_p->scl_xsts)
      {
         scale = evt_p->amp_sf;
      }
      else
      {
         scale = evt_p->event_status & 0x03; 
      }
      if (scale == 3)
      {
         scale++;
      }

      /* casting to perserve algorithm from proc_prd2epr code */ 
      if ((spha =(0.5 * (double)(scale * evt_p->sum_amps)/inp_p->amp_gain))!= 0)
      {
         double ratio = ((double)evt_p->pha) / ((double)spha);

         if ((ratio > (1 + inp_p->pha_ratio)) ||
             (ratio < (1 - inp_p->pha_ratio))) 
         {
            stat_p->bad_pha_ratio++; 
            evt_p->status |= HDET_PHA_RTO_STS; 
         } 
      }
      else
      {
         evt_p->status |= HDET_ZERO_PSUM_STS; 
      }
   }
}

