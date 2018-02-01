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


#include "hrc_process_events.h"

/*H***********************************************************************
 
* FILE NAME: hpe_setup_degap_file.c
 
* DESCRIPTION:
 
  The routine hpe_setup_degap_file() is called by hrc_process_events to
  populate the degap tables which will be used to transform the coarse 
  coordinates of the input event file to the chip coordinates of the 
  output event file.

* NOTES

  Errors are recorded by pushing them onto an error stack maintained 
  via the dserror lib. The calling routine can access the stack to
  determine if any errors occurred.
 
* REVISION HISTORY:
 
        Ref. No.        Date
        --------        ----
        1.0             20 Oct 1999
 
*H***********************************************************************/


void hpe_setup_degap_file(
   INPUT_PARMS_P_T inp_p,       /* I - input parameter data structure   */
   DEGAP_CONFIG_P_T*  dgp_p,    /* O - degap configuartion structure    */ 
   dsErrList*      hpe_err_p)   /* O - error stack pointer              */
{

   if (!dgp_p)
   {
      /* error- NULL pointer passed in */
      dsErrAdd(hpe_err_p, dsHBBNOMEMALLOCERR, Individual, Generic, 
               "dgp_p (degap configuration)");
   }
   else if (!(*dgp_p))
   {
      /* initialize degap table (NONE case is default config */
      l1h_configure_degap_table(inp_p->instrume, dgp_p, hpe_err_p); 
   
      if (*dgp_p)
      {
         /* copy data from input params structure to degap structure */
         (*dgp_p)->cf[0][0] = inp_p->cf[0][0];  
         (*dgp_p)->cf[0][1] = inp_p->cf[0][1];  
         (*dgp_p)->cf[1][0] = inp_p->cf[1][0];  
         (*dgp_p)->cf[1][1] = inp_p->cf[1][1];  
         strcpy((*dgp_p)->file, inp_p->degap_file); 

         /* overwrite NONE config case if necessary */
         if (!ds_strcmp_cis(inp_p->degap_file, "COEFF"))
         {
            l1h_set_coefficient_values((*dgp_p), 
                                       (*dgp_p)->cf[0][0], (*dgp_p)->cf[0][1],
                                       (*dgp_p)->cf[1][0], (*dgp_p)->cf[1][1],
                                       hpe_err_p);
         }
         else if (ds_strcmp_cis(inp_p->degap_file, "NONE")) 
         {
            /* load degapfile specified */
            l1h_load_degap_file(inp_p->degap_file, (*dgp_p), hpe_err_p);  
         } 
      } 
      else
      {
         /* error- allocation or setup failed */ 
         dsErrAdd(hpe_err_p, dsHPEDEGAPLOADERR, Individual, Generic);
      } 
   }
   else
   {
      /* warning degap table already set- not overriding */ 
      dsErrAdd(hpe_err_p, dsHBBDEGAPREREADERR, Individual, Generic);
   } 
}


