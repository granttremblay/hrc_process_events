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
 
* FILE NAME: t_hrc_process_events.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION:
 
  Unix main for the hrc_process_events program. This program is intended 
  to be an open-iraf C program to compute detector coordinates for hrc 
  event files. 

* NOTES:
 
  wmclaugh@cfa	March 25, 1996  First Version.
 
* REVISION HISTORY:
 
        Ref. No.        Date
        --------        ----
        1.1             25 Mar 1996
 
*H***********************************************************************/


#include <stdio.h>
#include <parameter.h>

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif


int main(int argc, char** argv)
{
  dsErrGroup groups_to_add_t = dsPTGRPERR;
  dsErrCode fail_status_t = dsNOERR;

  /* set the jump buffer, for the error lib.  The error library performs
     signal handling, and must have a place to jump to. */
  if(setjmp(dserr_jmpbuf) == 0)
  {

    /* init error handling routine */
    if((fail_status_t = dsErrInitLib(groups_to_add_t, argv[0])) == dsNOERR)
    {
      /* OPEN THE PRAMETER FILE */
      if(clinit(argv, argc, "rw") == NULL)
      {
         fail_status_t = dsOPENPARAMFERR;
         err_msg(dsOPENPARAMFSTDMSG, "hrc_process_events.par");
	 err_msg("ERROR: Parameter library error: %s.\n", paramerrstr());
      }
      else
      {    
         /* EXECUTE OUR PROGRAM */ 
         fail_status_t = hrc_process_events();
    
         /* CLOSE PARAMETER FILE AND RETURN TO THE OS */
         clclose();
      }

      dsErrCloseLib();
    }

  } /* end if(setjmp) */
  else
  {
    fail_status_t = dsGENERICERR;
  }


  return (fail_status_t); 
}
