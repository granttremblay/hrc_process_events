/*                                                                
**  Copyright (C) 1999-2007  Smithsonian Astrophysical Observatory 
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
*
* FILE NAME: hrc_process_configure_pixlib.c
*
* DEVELOPEMENT: tools
*
* JCC(2/2002) - pass geompar to the pixlib call.
* (6/2004)-condition check on pixlib
*H**************************************************************************/

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

void hrc_process_configure_pixlib(
   EVENT_SETUP_P_T   evtin_p,  /* I - input event file             */
   INPUT_PARMS_P_T   inp_p,    /* I/O - input parameters/data        */
   INST_KEYWORDS_P_T inst_p,   /* instrument keyword structure pntr */
   ALIGNMENT_REC_P_T aln_p,    /* I - ptr to an alignment file record */
   dsErrList*        err_p)    /* O   - error list pointer           */
{
   /* initialize pixlib for either flight or xrcf */
   if (inp_p->processing == HRC_PROC_FLIGHT)
   {
      if (pix_init_pixlib("FLIGHT",inp_p->geompar)!=PIX_GOOD)
      {
        /* fatal error */
         dsErrAdd(err_p, dsINITLIBERR, Individual, Custom,
         " ERROR: Could not initialize pixlib. Please check the setup.\n");
      }
   }
   else
   {
      if (pix_init_pixlib("XRCF",inp_p->geompar)!=PIX_GOOD)
      {
        /* fatal error */
         dsErrAdd(err_p, dsINITLIBERR, Individual, Custom,
         " ERROR: Could not initialize pixlib. Please check the setup.\n");
      }
   }
   if (err_p->contains_fatal== 0)
   { 
      inp_p->pix_init = TRUE;
      pix_set_randseed(inp_p->rand_seed);
 
      set_up_mirror(evtin_p->extension, inp_p,
                    evtin_p->file, aln_p, err_p);

      pix_get_fp_scale_in_asec(&inp_p->fp_scale);

      /* get default systems */
      pix_get_chipsys_default(inst_p->chip_system);
      pix_get_tdetsys_default(inst_p->tdet_system);
      pix_get_fpsys_default(inst_p->fp_system);
   }
}
