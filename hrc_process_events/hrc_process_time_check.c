/*                                                                
**  Copyright (C) 1998-2007  Smithsonian Astrophysical Observatory 
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


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif


/*************************************************************************

* DESCRIPTION

  The routine hrc_process_time_check() is called by hrc_process_events
  to determine whether the time range specified by TSTART and TSTOP in 
  the input obs.par file covers the time range of the input events. If 
  event times are outside of the range specified by the TSTART and TSTOP
  parameters an error warning is added to the error list. 

  The routine also prints the tstart ranges of the event files and the
  obs.par file if the debug level is set to 2 or higher. 

* NOTES 

  The routine does not return any error status, but detected errors are
  added to the error list passed in. 

*************************************************************************/

void hrc_process_time_check(
   INPUT_PARMS_P_T  inp_p,     /* I   - input parameters data structure */
   FILE*            log_p,     /* I   - pointer to debug log destination */
   dsErrList*       err_p)     /* O   - error list pointer              */
{

   /* only check if obs.par file used */
   if ((ds_strcmp_cis(inp_p->obsfile, "NONE") != 0) &&  
       (inp_p->obsfile[0] != '\0')) 
   {
      if (inp_p->evt_tstart < inp_p->obs_tstart)   
      {
         dsErrAdd(err_p, dsAPETSTARTERR, Individual, Generic, 
                  inp_p->evt_tstart, inp_p->obsfile, inp_p->obs_tstart);
      } 
      if (inp_p->evt_tstop > inp_p->obs_tstop)   
      {
         dsErrAdd(err_p, dsAPETSTOPERR, Individual, Generic, 
                  inp_p->evt_tstop, inp_p->obsfile, inp_p->obs_tstop);
      } 
      if (inp_p->obs_tstop < inp_p->obs_tstart)
      {
         dsErrAdd(err_p, dsAPEOBSTIMERNGERR, Individual, Generic, 
                  inp_p->obs_tstop, inp_p->obs_tstart, inp_p->obsfile);
      } 

      if (inp_p->debug > DEBUG_LEVEL_1)
      {
         fprintf(log_p, "\n Time interval specified in %s \n", 
            inp_p->obsfile); 
         fprintf(log_p, "    TSTART  : %12.6f\n", inp_p->obs_tstart); 
         fprintf(log_p, "    TSTOP   : %12.6f\n", inp_p->obs_tstop); 
      } 
   } 

   if ((inp_p->debug > DEBUG_LEVEL_1) && (inp_p->evt_tstart != DBL_MAX) &&
       (inp_p->evt_tstop != DBL_MIN))
   {
      fprintf(log_p, "\n Time interval of events in event file(s)\n"); 
      fprintf(log_p, "    TSTART  : %12.6f\n", inp_p->evt_tstart); 
      fprintf(log_p, "    TSTOP   : %12.6f\n", inp_p->evt_tstop); 
   } 
}  

