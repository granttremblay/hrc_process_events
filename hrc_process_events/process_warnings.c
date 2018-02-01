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

/*H***********************************************************************
 
* FILE NAME: process_warnings.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION:

  The routine process_warnings is called by hrc_process_events to write
  out any warning messages which may have been placed onto the error queue
  during the processing of an input file. The routine returns a boolean
  value which is set to TRUE if the list contains a warning message which 
  would cause the input file to be considered a 'bad' file for sttistical 
  purposes. The routine removes the warning messages from the input error
  list but does nothing to the fatal errors.

  NOTE: At present there are no warning messages which would cause the 
  input file to be considered bad (these instances are all fatal errors).
  The code has been left in place for future use.  
 
*H***********************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 


boolean process_warnings(dsErrList* all_err_p, FILE* log_ptr, int level)
{
   dsErrList* warnings_p = NULL;
   boolean increment_badfile = FALSE; 

   dsErrCreateList(&warnings_p);

   /* copy warnings from main error list to warning list */
   dsErrLookAllSev(all_err_p, Warning, warnings_p);
 
   /* remove warnings from main error list */
   dsErrRemoveAllSev(all_err_p, Warning);

   /* print out warning messages */
   dsErrPrintList(warnings_p, dsErrTrue);

   /* write warnings to logfile, if specified */
   if ((level > DEBUG_LEVEL_0) && (log_ptr != stderr) && (log_ptr != stdout))
   {
      dsErrDirectOutput(log_ptr);
      dsErrPrintList(warnings_p, dsErrTrue);
      dsErrDirectOutput(stderr);
   }

   /* search for warnings which make event file 'bad'
    * and set flag to TRUE if found- currently no such
    * warnings are defined 
    */ 

   /* remove the warning messages */
   dsErrRemoveAllSev(warnings_p, Warning);
 
   return (increment_badfile);  

}





