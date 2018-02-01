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
 
* FILE NAME: hrc_process_configure_pixlib.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: 
 
* NOTES:

*H**************************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 


void hrc_process_set_instrume(
   INPUT_PARMS_P_T   inp_p,    /* I/O - input parameters/data        */
   dsErrList*        err_p)    /* O   - error list pointer           */
{
   int err = 0; 
   int pos;   /* position in instrume character string */
   char   default_instrume[] = "hrc-s"; 
   char*  instrument = default_instrume; 

   /* need to convert instrume(nt) keyword to lowercase */
   if ((pos = strlen(inp_p->instrume)) > DS_SZ_KEYWORD)
   {
      pos = DS_SZ_KEYWORD;
   }
   while (pos-- > 0)
   {
      inp_p->instrume[pos] = (char) tolower(inp_p->instrume[pos]);
   }

   if (strcmp(inp_p->instrume,"hrc-i") == 0)
   {
      inp_p->hrc_system = HRC_IMG_SYS;
      instrument = inp_p->instrume; 
   }
   else if (strcmp(inp_p->instrume,"hrc-s") == 0)
   {
      inp_p->hrc_system = HRC_SPC_SYS;
      instrument = inp_p->instrume; 
   }
   else if (strcmp(inp_p->instrume,"hsi") == 0)
   {
      inp_p->hrc_system = HSI_IMG_SYS;
      instrument = inp_p->instrume; 
   }
   else if (strcmp(inp_p->instrume,"hrc-si") == 0)
   {
      inp_p->hrc_system = HRC_SPC_IMG_SYS;
   }
   else
   {
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Custom,
               "ERROR: The detector specified by the obs.par detnam or the hrc_process_events.par parameter instrume is not supported."); 
      err = 1; 
   }
}
