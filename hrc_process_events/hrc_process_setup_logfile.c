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
 
* DESCRIPTION:  The routine hrc_process_setup_logfile is called to 
  initialize the output debug log file. This routine uses the logfile 
  parameter stored in INPUT_PARMS_T to determine whether to echo debugging
  statements to stdout or to an external file. If the parameter is set to 
  "stdout" the comments are output to stdout (typically the terminal 
  screen).
 
  The routine returns a file descriptor which is either set to an external
  file or to stdout. If a value of "NONE" or "blank" is specified no debug
  info is written out (in effect it overwrites the debug value).
 
  NOTE:  The debug value as specified in the parameter file and stored in
  INPUT_PARMS_T determines the amount of information displayed. A value of
  0 indicates no information. A value of 5 indicates detailed logging
 
*************************************************************************/
 
 
FILE* hrc_process_setup_logfile(
   INPUT_PARMS_P_T inp_p, /* I - input parameters */
   dsErrList*      err_p) /* O - pointer to error message stack */
{
   FILE *log_p = NULL;
 
   if ((inp_p->debug > DEBUG_LEVEL_0)&&
        (ds_strcmp_cis(inp_p->logfile, HDET_LOG_TO_STDOUT) != 0))
    {
       if ((ds_strcmp_cis(inp_p->logfile, "none") == 0) ||
           (strcmp(inp_p->logfile, "\0") == 0))
       {
          inp_p->debug = 0; 
       } 
       else if ((log_p = fopen(inp_p->logfile, "w")) == NULL)
       {
          /* if can't open log file- sent debug info to stdout */
          log_p = stdout;
          dsErrAdd(err_p, dsHPELOGOPENERR, Individual, Generic,
                   inp_p->logfile);
       }
    }
    else
    {
       log_p = stdout;
    }

    if (inp_p->debug > DEBUG_LEVEL_0)
    {
       fprintf(log_p, " input stack    : %s\n", inp_p->stack_in);
       fprintf(log_p, " output qp file : %s\n", inp_p->outfile);
       fprintf(log_p, " alignment file : %s\n", inp_p->align_file);
       fprintf(log_p, " aspect file    : %s\n", inp_p->asp_file);
       fprintf(log_p, " adc file       : %s\n", inp_p->adc_file);
       fprintf(log_p, " gain file      : %s\n", inp_p->gain_file);
       fprintf(log_p, " degap file     : %s\n", inp_p->degap_file);
       fprintf(log_p, " observation file : %s\n", inp_p->obsfile);
       fprintf(log_p, " bad pixel file : %s\n", inp_p->badpixfile);
    } 
    if (inp_p->debug > DEBUG_LEVEL_1)
    {
       fprintf(log_p, " output columns : %s\n", inp_p->outcols);
       fprintf(log_p, " grid ratio     : %7.4f\n", inp_p->grid_ratio);
       fprintf(log_p, " pha ratio      : %7.4f\n", inp_p->pha_ratio);
       fprintf(log_p, " wire charge    : %4d\n", inp_p->wire_charge);
       fprintf(log_p, " x corr factor 1: %7.4f\n",
               inp_p->cf[HDET_PLANE_X][HDET_1ST_ORD_CF] );
       fprintf(log_p, " x corr factor 2: %7.4f\n",
               inp_p->cf[HDET_PLANE_X][HDET_2ND_ORD_CF] );
       fprintf(log_p, " y corr factor 1: %7.4f\n",
               inp_p->cf[HDET_PLANE_Y][HDET_1ST_ORD_CF] );
       fprintf(log_p, " y corr factor 2: %7.4f\n",
               inp_p->cf[HDET_PLANE_Y][HDET_2ND_ORD_CF] );
       fprintf(log_p, " debug level    : %1d\n\n", inp_p->debug);
    }
 
   return(log_p);
}
 
 
