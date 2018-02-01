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
 
* FILE NAME: write_instrume_params.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION:
 
  The file write_instrume_params.c contains the following modules used by
  hrc_process_events() in dealing with instrument paramater keywords:
 
        read_instrume_params() 
        write_instrume_params() 
 
* NOTES:
 
  wmclaugh@cfa  December 06, 1996  First Version.
 
* REVISION HISTORY:
 
        Ref. No.        Date
        --------        ----
        1.1             06 dec 1996
*H***********************************************************************/


#ifndef HRC_PROCESS_EVENTS_H 
#include "hrc_process_events.h"
#endif

/*************************************************************************
 
* DESCRIPTION
 
  The routine write_instrume_params() writes parameters from the 
  instrument keyword data structure to the specified output event file. 
 
*************************************************************************/

void write_instrume_params(
   dmBlock*          extension,   /* I - output file name */
   INST_KEYWORDS_P_T inst_p)
{
    char key_val[DS_SZ_KEYWORD];

    sprintf(key_val, "RAW:%s", inst_p->chip_system ); 
    dmKeyWrite_c(extension, "ACSYS1", key_val, NULL,
                 "reference for raw chip coord system");
    sprintf(key_val, "CHIP:%s", inst_p->chip_system ); 
    dmKeyWrite_c(extension, "ACSYS2", key_val, NULL,
                 "reference for degap corrected chip coord system");
    sprintf(key_val, "TDET:%s", inst_p->tdet_system ); 
    dmKeyWrite_c(extension, "ACSYS3", key_val, NULL,
                 "reference for tiled detector coord system");
    sprintf(key_val, "DET:%s", inst_p->fp_system ); 
    dmKeyWrite_c(extension, "ACSYS4", key_val, NULL,
                 "reference for focal plane coord system");
    sprintf(key_val, "SKY:%s", inst_p->fp_system ); 
    dmKeyWrite_c(extension, "ACSYS5", key_val, NULL,
                 "reference for sky coord system");
}


