/*                                                                
**  Copyright (C) 1998-2010  Smithsonian Astrophysical Observatory 
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

/*************************************************************************

* FILE NAME: hrc_setup_columns.c
 
* DESCRIPTION:

  The routine hrc_setup_columns() is called upon by hrc_process_events to
  create the columns of the output event file. It is responsible for 
  obtaining the column descriptor as well as for setting the values of
  TLMIN/TLMAX/TUNITS keywords as specified in the level 1 hrc to archive
  icd. (Macros defining these values are found in ds_hrc_config.h).

* NOTES

  The routine does not pass a return status back, but detected errors are
  added to the maintained error list which is passed in as a parameter.  


 JCC(7/8/02) - change (theta,phi) to (phi,theta)
 JCC(11/18/02)-fix 'ciao2.3 h_p_e issue' { incorrect reference points 
               of sky_XY (EQPOS[RA,DEC])}
    ( when obsfile=NONE, get crval[0..1] from ra_nom/dec_nom of EVT1.) 
    ( fix: first, set ra_nom/dec_nom properly in hrc_process_read_obsfile 
      and hrc_process_setup_input_file, then assign them to crval[0..1] 
      in hrc_setup_columns. )
 JCC(2/2003)- remove EventKey that won't be defined in new hdrlib.
            - remove "event_hdr.h"
1/2009-for old/new gain img, HDET_PI: maxVal=255/512, dtype=short/long (obsolete10/2009)
10/2009- replace peter hrcS 3dim gain image w/ dph new hrcS gain table.
10/2009- add fap new hrcI gain image.
       - see 'Notes on outCol PI'
1/2010 - fix dong-woo's ontime/GTI_CPT2 problem which was actually caused by 
         the subspace bugs in the datamodel library. 
         Add hpeSetRange_s and hpe_set_ranges as a workaround.
5/2010 - fix pi subspace problem w/ new gainfile
*************************************************************************/

#include <ctype.h>
#include <limits.h>

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#define HRC_PROCESS_EVENTS_H
#endif

#ifndef DSLIB_H
#include "dslib.h"
#define DSLIB_H
#endif

#ifndef DS_HRC_CONFIG_H
#include "ds_hrc_config.h"
#define DS_HRC_CONFIG_H
#endif

void hrc_setup_columns(
   EVENT_SETUP_P_T evtout_p,      /* I/O - output file inprocession     */
   INPUT_PARMS_P_T inp_p,         /* I/O - input parameters/data        */
   char**          out_names,     /* I   - output event column names    */ 
   dsErrList*      hpe_err_p)     /* O   - error list pointer           */

{ 
   short rr;

   double sky_x_min=0.0, sky_x_max=0.0, sky_y_min=0.0, sky_y_max=0.0;
   double det_x_min=0.0, det_x_max=0.0, det_y_min=0.0, det_y_max=0.0;



   for (rr = 0; rr < evtout_p->num_cols; rr++)
   {
      switch (evtout_p->mapping[rr])
      {
         case HDET_SKY:    /* vectored column pair- X,Y */  
         {
            VEC2_LONG axlen;
            char      names[2][DS_SZ_KEYWORD];
            char*     cptnames[2];
            char*     eqnames[] = {"RA", "DEC"};
            dmDescriptor* temp_desc;
   
            cptnames[0] = names[0]; 
            cptnames[1] = names[1]; 
            strcpy(cptnames[0], "X"); 
            strcpy(cptnames[1], "Y"); 
            /* if output column name is lowercase, set x/y to lc */
            if (islower((int)out_names[rr][0]))
            {
               string_to_lowercase(cptnames[0], DS_SZ_KEYWORD);
               string_to_lowercase(cptnames[1], DS_SZ_KEYWORD);
            }
   
            pix_get_fp_image_cntr(inp_p->crpix);

            /* 11/2002 */
            /* EventKey*  obs_info_p  - ra/dec_nom from observation.par */
            /* inp_p->crval[0] = obs_info_p->ObservKey_ptr.ra_nom; */
            /* inp_p->crval[1] = obs_info_p->ObservKey_ptr.dec_nom; */
            inp_p->crval[0] = inp_p->ra_nom;
            inp_p->crval[1] = inp_p->dec_nom;

            inp_p->crval[2] = 0.0;    

            inp_p->cdelt[0] = inp_p->fp_scale/-3600.0;
            inp_p->cdelt[1] = inp_p->cdelt[0]*-1;

            evtout_p->dim[rr] = 2;
            evtout_p->desc[rr] =
               dmColumnCreateVector(evtout_p->extension, out_names[rr],
                                    evtout_p->types[rr], 0, 
                                    DS_HRC_SKY_TUNITS, "Sky coords",
                                    cptnames, evtout_p->dim[rr]);
            dmCoordCreate_d(evtout_p->desc[rr], "EQPOS", "deg",
                            eqnames, 2, "TAN", inp_p->crpix, inp_p->crval, 
                            inp_p->cdelt, NULL);
    
            /* set TLMIN/TLMAX for sky coords */
            pix_get_fp_image_dimn(axlen);

	    switch (inp_p->hrc_system)
               {
	       case HRC_IMG_SYS:
		 sky_x_min=DS_HRCI_X_TLMIN; sky_x_max=DS_HRCI_X_TLMAX;
		 sky_y_min=DS_HRCI_Y_TLMIN; sky_y_max=DS_HRCI_Y_TLMAX;
		 break;
	       case HRC_SPC_SYS:
	       case HRC_SPC_IMG_SYS:
		 sky_x_min=DS_HRCS_X_TLMIN; sky_x_max=DS_HRCS_X_TLMAX;
		 sky_y_min=DS_HRCS_Y_TLMIN; sky_y_max=DS_HRCS_X_TLMAX;
		 break;
	       case HSI_IMG_SYS:
		 sky_x_min=DS_HRCI_X_TLMIN; sky_x_max=axlen[0];
		 sky_y_min=DS_HRCI_Y_TLMIN; sky_y_max=axlen[1];
		 break;

	       default:
		 break;
	       }
            switch (evtout_p->types[rr])
            {

               case dmLONG:
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_l(temp_desc, 
                                            sky_x_min, sky_x_max);
                  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_l(temp_desc, 
                                            sky_y_min, sky_y_max);
                  }
	       

               break;
               case dmSHORT:
               {
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_s(temp_desc, 
					    sky_x_min, sky_x_max);
                  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_s(temp_desc, 
					    sky_y_min, sky_y_max);
                  }
               }
               break;
               case dmFLOAT:
               {
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_f(temp_desc, 
					    sky_x_min, sky_x_max);
		  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_f(temp_desc, 
					    sky_y_min, sky_y_max );
                  }
               }
               break;
               case dmDOUBLE:
               {
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_d(temp_desc,
					    sky_x_min, sky_x_max );
                  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_d(temp_desc, 
					    sky_y_min, sky_y_max );
		  }
               }
               break;
               default:
                  dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         }
         break; 
         case HDET_RAW:   /* vectored coord pair- RAWX/RAWY */ 
         {
            char  names[2][DS_SZ_KEYWORD];
            char* cptnames[2];
            dmDescriptor* temp_desc;

            cptnames[0] = names[0]; 
            cptnames[1] = names[1]; 
            strcpy(cptnames[0], HDET_RAW_X_COL); 
            strcpy(cptnames[1], HDET_RAW_Y_COL); 

            /* if output column name is lowercase, set x/y to lc */
            if (islower((int)out_names[rr][0]))
            {
               string_to_lowercase(cptnames[0], DS_SZ_KEYWORD);
               string_to_lowercase(cptnames[1], DS_SZ_KEYWORD);
            }

            inp_p->do_raw = TRUE; 
            evtout_p->dim[rr] = 2;
            evtout_p->desc[rr] =
               dmColumnCreateVector(evtout_p->extension, out_names[rr],
                                    evtout_p->types[rr], 0,
                                    DS_HRC_RAW_TUNITS, "Raw coords",
                                    cptnames, evtout_p->dim[rr]);

            if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))  
            {
               switch (inp_p->hrc_system)
               {
                  case HRC_IMG_SYS:
                     dmDescriptorSetRange_l(temp_desc,
                        DS_HRCI_RAWX_TLMIN, DS_HRCI_RAWX_TLMAX);
                  break;
                  case HRC_SPC_SYS:
                     dmDescriptorSetRange_l(temp_desc,
                        DS_HRCS_RAWX_TLMIN,DS_HRCS_RAWX_TLMAX);
                  break;
                  case HRC_SPC_IMG_SYS:
                     dmDescriptorSetRange_l(temp_desc,
                        DS_HRCS_RAWX_TLMIN,DS_HRCS_RAWX_TLMAX);
                  break;
                  case HSI_IMG_SYS:
                     dmDescriptorSetRange_l(temp_desc,
                        DS_HRCI_RAWX_TLMIN, DS_HRCI_RAWX_TLMAX);
                  break;
                  default:
                     /* error - unknown system type */
                     dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))  
            {
               switch (inp_p->hrc_system)
               {
                  case HRC_IMG_SYS:
                     dmDescriptorSetRange_l(temp_desc,
                        DS_HRCI_RAWY_TLMIN, DS_HRCI_RAWY_TLMAX);
                  break;
                  case HRC_SPC_SYS:
                     dmDescriptorSetRange_l(temp_desc,
                        DS_HRCS_RAWY_TLMIN,DS_HRCS_RAWY_TLMAX);
                  break;
                  case HRC_SPC_IMG_SYS:
                     dmDescriptorSetRange_l(temp_desc, 
                        DS_HRCS_RAWY_TLMIN,DS_HRCS_RAWY_TLMAX);
                  break;
                  case HSI_IMG_SYS:
                     dmDescriptorSetRange_l(temp_desc,
                        DS_HRCI_RAWY_TLMIN, DS_HRCI_RAWY_TLMAX);
                  break;
                  default:
                     /* error - unknown system type */
                     dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
         } 
         break; 
         case HDET_CHIP:   /* vectored coord pair- CHIPX/CHIPY */ 
         {
            double crpix[2] = {0.5, 0.5},
                   crval[2] = {0.0, 0.0},
                   cdelt[2];
            char   names[2][DS_SZ_KEYWORD];
            char*  cpcnames[] = {"CPCX", "CPCY"};
            char*  cptnames[2];
            dmDescriptor* temp_desc;

            cptnames[0] = names[0]; 
            cptnames[1] = names[1]; 
            strcpy(cptnames[0], HDET_CHIP_X_COL); 
            strcpy(cptnames[1], HDET_CHIP_Y_COL); 
 
            /* if output column name is lowercase, set x/y to lc */
            if (islower((int)out_names[rr][0]))
            {
               string_to_lowercase(cptnames[0], DS_SZ_KEYWORD);
               string_to_lowercase(cptnames[1], DS_SZ_KEYWORD);
            }

            evtout_p->dim[rr] = 2;
            evtout_p->desc[rr] =
               dmColumnCreateVector(evtout_p->extension, out_names[rr],
                                    evtout_p->types[rr], 0,
                                    DS_HRC_CHIP_TUNITS, "Chip coords",
                                    cptnames, evtout_p->dim[rr]);

            pix_get_chip_scale_in_mm(cdelt);
            dmCoordCreate_d(evtout_p->desc[rr], "CPC", "mm",
                            cpcnames, 2, "LINEAR", crpix, crval, cdelt,
                            NULL);

            if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
            {
               switch (inp_p->hrc_system)
               {
                  case HRC_IMG_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCI_CHIPX_TLMIN, DS_HRCI_CHIPX_TLMAX);
                  break;
                  case HRC_SPC_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCS_CHIPX_TLMIN,DS_HRCS_CHIPX_TLMAX);
                  break;
                  case HRC_SPC_IMG_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCS_CHIPX_TLMIN,DS_HRCS_CHIPX_TLMAX);
                  break;
                  case HSI_IMG_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCI_CHIPX_TLMIN, DS_HRCI_CHIPX_TLMAX);
                  break;
                  default:
                     /* error - unknown system type */
                     dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
            {
               switch (inp_p->hrc_system)
               {
                  case HRC_IMG_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCI_CHIPY_TLMIN, DS_HRCI_CHIPY_TLMAX);
                  break;
                  case HRC_SPC_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCS_CHIPY_TLMIN,DS_HRCS_CHIPY_TLMAX);
                  break;
                  case HRC_SPC_IMG_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCS_CHIPY_TLMIN,DS_HRCS_CHIPY_TLMAX);
                  break;
                  case HSI_IMG_SYS:
                     dmDescriptorSetRange_s(temp_desc,
                        DS_HRCI_CHIPY_TLMIN, DS_HRCI_CHIPY_TLMAX);
                  break;
                  default:
                     /* error - unknown system type */
                     dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
         }
         break; 
         case HDET_TDET:   /* vectored coord pair- TDETX/TDETY */ 
         {
            char  names[2][DS_SZ_KEYWORD];
            char* cptnames[2];
            dmDescriptor* temp_desc;

            cptnames[0] = names[0]; 
            cptnames[1] = names[1]; 
            strcpy(cptnames[0], HDET_TDET_X_COL); 
            strcpy(cptnames[1], HDET_TDET_Y_COL); 

            /* if output column name is lowercase, set x/y to lc */
            if (islower((int)out_names[rr][0]))
            {
               string_to_lowercase(cptnames[0], DS_SZ_KEYWORD);
               string_to_lowercase(cptnames[1], DS_SZ_KEYWORD);
            }

            evtout_p->dim[rr] = 2;
            evtout_p->desc[rr] =
               dmColumnCreateVector(evtout_p->extension, out_names[rr],
                                    evtout_p->types[rr], 0,
                                    DS_HRC_TDET_TUNITS, "Tdet coords",
                                    cptnames, evtout_p->dim[rr]);
            if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
            {
               double tlmax_value;
               VEC2_LONG tdet_dim;

               if (pix_get_tdet_image_dimn(tdet_dim) == PIX_GOOD)
               {
                  tlmax_value = tdet_dim[0]; 
               }
               else if ((inp_p->hrc_system == HRC_IMG_SYS) || 
                   (inp_p->hrc_system == HSI_IMG_SYS))
               {
                  tlmax_value = DS_HRC_TDETX_TLMAX;  
               }
               else 
               {
                  tlmax_value = DS_HRCS_TDETX_TLMAX;  
               }

               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(temp_desc,
                        DS_HRC_TDETX_TLMIN, tlmax_value);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(temp_desc,
                        (float) DS_HRC_TDETX_TLMIN, 
                        (float) tlmax_value);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(temp_desc,
                        (long) DS_HRC_TDETX_TLMIN, (long) tlmax_value);
                  break;
                  case dmSHORT:
                  {
                     short cast_min = (short) DS_HRC_TDETX_TLMIN,
                           cast_max = (short) tlmax_value;
    
                     if (DS_HRC_TDETX_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRC_TDETX_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(temp_desc,
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            } 
            if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
            {
               double tlmax_value;
               VEC2_LONG tdet_dim;

               if (pix_get_tdet_image_dimn(tdet_dim) == PIX_GOOD)
               {
                  tlmax_value = tdet_dim[1]; 
               }
               else if ((inp_p->hrc_system == HRC_IMG_SYS) ||
                   (inp_p->hrc_system == HSI_IMG_SYS))
               {
                  tlmax_value = DS_HRC_TDETY_TLMAX;
               }
               else
               {
                  tlmax_value = DS_HRCS_TDETY_TLMAX;
               }

               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(temp_desc,
                        DS_HRC_TDETY_TLMIN, tlmax_value);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(temp_desc,
                        (float) DS_HRC_TDETY_TLMIN, (float) tlmax_value);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(temp_desc,
                        (long) DS_HRC_TDETY_TLMIN, (long) tlmax_value);
                  break;
                  case dmSHORT:
                  {
                     short cast_min = (short) DS_HRC_TDETY_TLMIN,
                           cast_max = (short) tlmax_value;
   
                     if (DS_HRC_TDETY_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRC_TDETY_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(temp_desc,
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
         }
         break;
         case HDET_DET:   /* vectored coord pair- DETX/DETY */ 
         {
            double    crpix[2],
                      crval[2] = {0.0, 0.0},
                      cdelt[2];
            double    tanp_params[4] = {0.0,     /* Rotation (deg)         */
                                        270.0,   /* pole longitude (deg)   */
                                        90.0,    /* pole latitude (deg)    */
                                        2000.0}; /* equinox (Julian years) */
            VEC2_LONG axlen;
            char      names[2][DS_SZ_KEYWORD];
            char*     cptnames[2];
            char*     mscnames[] = {"PHI", "THETA"};  /* 7/8/02 */
            dmDescriptor* temp_desc;
 
            cptnames[0] = names[0]; 
            cptnames[1] = names[1]; 
            strcpy(cptnames[0], HDET_DET_X_COL); 
            strcpy(cptnames[1], HDET_DET_Y_COL); 

            /* if output column name is lowercase, set x/y to lc */
            if (islower((int)out_names[rr][0]))
            {
               string_to_lowercase(cptnames[0], DS_SZ_KEYWORD);
               string_to_lowercase(cptnames[1], DS_SZ_KEYWORD);
            }

            evtout_p->dim[rr] = 2;
            evtout_p->desc[rr] =
               dmColumnCreateVector(evtout_p->extension, out_names[rr],
                                    evtout_p->types[rr], 0,
                                    DS_HRC_DET_TUNITS, "Det coords",
                                    cptnames, evtout_p->dim[rr]);

            pix_get_fp_image_cntr(crpix);

            cdelt[0] = cdelt[1] = inp_p->fp_scale/3600.0;
            dmCoordCreate_d(evtout_p->desc[rr], "MSC", "deg",
                            mscnames, 2, "TAN-P", crpix, crval, cdelt,
                            tanp_params);

            pix_get_fp_image_dimn(axlen);
	    switch (inp_p->hrc_system)
               {
	       case HRC_IMG_SYS:
		 det_x_min=DS_HRCI_X_TLMIN; det_x_max=DS_HRCI_X_TLMAX;
		 det_y_min=DS_HRCI_Y_TLMIN; det_y_max=DS_HRCI_Y_TLMAX;
		 break;
	       case HRC_SPC_SYS:
	       case HRC_SPC_IMG_SYS:
		 det_x_min=DS_HRCS_X_TLMIN; det_x_max=DS_HRCS_X_TLMAX;
		 det_y_min=DS_HRCS_Y_TLMIN; det_y_max=DS_HRCS_X_TLMAX;
		 break;
	       case HSI_IMG_SYS:
		 det_x_min=DS_HRCI_X_TLMIN; det_x_max=axlen[0];
		 det_y_min=DS_HRCI_Y_TLMIN; det_y_max=axlen[1];
		 break;

	       default:
		 break;
	       }

            switch (evtout_p->types[rr])
            {
               case dmLONG:
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_l(temp_desc,
                                            det_x_min, det_x_max );
                  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_l(temp_desc,
                                            det_y_min, det_y_max );
                  }
               break;
               case dmSHORT:
               {
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_s(temp_desc, 
					    det_x_min, det_x_max );
                  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_s(temp_desc, 
					    det_y_min, det_y_max );
                  }
               }
               break;
               case dmFLOAT:
               {
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_f(temp_desc, 
					    det_x_min, det_x_max );
                  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_f(temp_desc,
					    det_y_min, det_y_max );
                  }
               }
               break;
               case dmDOUBLE:
               {
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 1)))
                  {
                     dmDescriptorSetRange_d(temp_desc,
					    det_x_min, det_x_max );
                  }
                  if ((temp_desc = dmGetCpt(evtout_p->desc[rr], 2)))
                  {
                     dmDescriptorSetRange_d(temp_desc,
					    det_y_min, det_y_max );
                  }
               }
               break;
               default:
                  dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         } 
         break; 
         case HDET_TIME:
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                              evtout_p->types[rr],  0, 
                                DS_HRC_TIME_TUNITS, "Time tag (TT)");
         break; 

         case HDET_PHA:      
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                              evtout_p->types[rr],  0, 
                              DS_HRC_PHA_TUNITS, 
                                "Pulse height"); 
         break; 

        /* -------------- */
         case HDET_PI:

          /*11/2009-outCol PI: maxPI=255||1023, SHORT dtype for old+new gainmap*/
          /*if ( inp_p->gainflag==OLD_SI_GAIN)   */
            evtout_p->types[rr] = dmSHORT;/*outCol PI: maxPI=255||1023*/

            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                              evtout_p->types[rr],  0, 
                              DS_HRC_PI_TUNITS, "Pulse Invariant"); 
         break; 
        /* -------------- */


         case HDET_RAW_X:    
            inp_p->do_raw = TRUE; 
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                              evtout_p->types[rr],  0, 
                              DS_HRC_RAW_TUNITS, "Raw X (no degapping)");
            switch (inp_p->hrc_system)
            {
               case HRC_IMG_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCI_RAWX_TLMIN, DS_HRCI_RAWX_TLMAX);
               break;
               case HRC_SPC_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCS_RAWX_TLMIN,DS_HRCS_RAWX_TLMAX);
               break;
               case HRC_SPC_IMG_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCS_RAWX_TLMIN,DS_HRCS_RAWX_TLMAX);
               break;
               case HSI_IMG_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCI_RAWX_TLMIN, DS_HRCI_RAWX_TLMAX);
               break;
               default:
                  /* error - unknown system type */
                  dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         break; 
         case HDET_RAW_Y:    
            inp_p->do_raw = TRUE; 
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, 
                             DS_HRC_RAW_TUNITS, "Raw Y (no degapping)");
            switch (inp_p->hrc_system)
            {
               case HRC_IMG_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCI_RAWY_TLMIN, DS_HRCI_RAWY_TLMAX);
               break;
               case HRC_SPC_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCS_RAWY_TLMIN,DS_HRCS_RAWY_TLMAX);
               break;
               case HRC_SPC_IMG_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCS_RAWY_TLMIN,DS_HRCS_RAWY_TLMAX);
               break;
               case HSI_IMG_SYS:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     DS_HRCI_RAWY_TLMIN, DS_HRCI_RAWY_TLMAX);
               break;
               default:
                  /* error - unknown system type */
                  dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         break; 
         case HDET_CHIP_X:   
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                              evtout_p->types[rr],  0,
                              DS_HRC_CHIP_TUNITS, 
                              "Chip X (degapping applied)");
            switch (inp_p->hrc_system)
            {
               case HRC_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCI_CHIPX_TLMIN, DS_HRCI_CHIPX_TLMAX);
               break;
               case HRC_SPC_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCS_CHIPX_TLMIN,DS_HRCS_CHIPX_TLMAX);
               break;
               case HRC_SPC_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCS_CHIPX_TLMIN,DS_HRCS_CHIPX_TLMAX);
               break;
               case HSI_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCI_CHIPX_TLMIN, DS_HRCI_CHIPX_TLMAX);
               break;
               default:
                  /* error - unknown system type */
                  dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         break;
         case HDET_CHIP_Y:   
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0,
                             DS_HRC_CHIP_TUNITS, 
                             "Chip Y (degapping applied)");
            switch (inp_p->hrc_system)
            {
               case HRC_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCI_CHIPY_TLMIN, DS_HRCI_CHIPY_TLMAX);
               break;
               case HRC_SPC_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCS_CHIPY_TLMIN,DS_HRCS_CHIPY_TLMAX);
               break;
               case HRC_SPC_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCS_CHIPY_TLMIN,DS_HRCS_CHIPY_TLMAX);
               break;
               case HSI_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCI_CHIPY_TLMIN, DS_HRCI_CHIPY_TLMAX);
               break;
               default:
                  /* error - unknown system type */
                  dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         break;
         case HDET_TDET_X:  
         {
            double tlmax_value;
            VEC2_LONG tdet_dim;

            if (pix_get_tdet_image_dimn(tdet_dim) == PIX_GOOD)
            {
               tlmax_value = tdet_dim[0]; 
            }
            else if ((inp_p->hrc_system == HRC_IMG_SYS) ||
                (inp_p->hrc_system == HSI_IMG_SYS))
            {
               tlmax_value = DS_HRC_TDETX_TLMAX;
            }
            else
            {
               tlmax_value = DS_HRCS_TDETX_TLMAX;
            }

            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0,
                             DS_HRC_TDET_TUNITS, "Tiled X");
            switch(evtout_p->types[rr])
            {
               case dmDOUBLE:
                  dmDescriptorSetRange_d(evtout_p->desc[rr],
                     DS_HRC_TDETX_TLMIN, tlmax_value);
               break;
               case dmFLOAT:
                  dmDescriptorSetRange_f(evtout_p->desc[rr],
                     (float) DS_HRC_TDETX_TLMIN, (float) tlmax_value);
               break; 
               case dmLONG:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     (long) DS_HRC_TDETX_TLMIN, (long) tlmax_value);
               break;
               case dmSHORT: 
               {
                  short cast_min = (short) DS_HRC_TDETX_TLMIN,
                        cast_max = (short) tlmax_value;
 
                  if (DS_HRC_TDETX_TLMIN < SHRT_MIN)
                  {
                     /* error - underflow */
                     cast_min = SHRT_MIN;
                  }
                  if (DS_HRC_TDETX_TLMAX > SHRT_MAX)
                  {
                     /* error - overflow */
                     cast_max = SHRT_MAX;
                  }
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     cast_min, cast_max);
               }
               break;
               default:
                  /* error invalid type */
                  dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         }
         break;
         case HDET_TDET_Y:  
         {
            double tlmax_value;
            VEC2_LONG tdet_dim;

            if (pix_get_tdet_image_dimn(tdet_dim) == PIX_GOOD)
            {
               tlmax_value = tdet_dim[1]; 
            }
            else if ((inp_p->hrc_system == HRC_IMG_SYS) ||
                (inp_p->hrc_system == HSI_IMG_SYS))
            {
               tlmax_value = DS_HRC_TDETY_TLMAX;
            }
            else
            {
               tlmax_value = DS_HRCS_TDETY_TLMAX;
            }

            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0,
                             DS_HRC_TDET_TUNITS, "Tiled Y");
            switch(evtout_p->types[rr])
            {
               case dmDOUBLE:
                  dmDescriptorSetRange_d(evtout_p->desc[rr],
                     DS_HRC_TDETY_TLMIN, tlmax_value);
               break;
               case dmFLOAT:
                  dmDescriptorSetRange_f(evtout_p->desc[rr],
                     (float) DS_HRC_TDETY_TLMIN, (float) tlmax_value);
               break; 
               case dmLONG:
                  dmDescriptorSetRange_l(evtout_p->desc[rr],
                     (long) DS_HRC_TDETY_TLMIN, (long) tlmax_value);
               break;
               case dmSHORT: 
               {
                  short cast_min = (short) DS_HRC_TDETY_TLMIN,
                        cast_max = (short) tlmax_value;
 
                  if (DS_HRC_TDETY_TLMIN < SHRT_MIN)
                  {
                     /* error - underflow */
                     cast_min = SHRT_MIN;
                  }
                  if (DS_HRC_TDETY_TLMAX > SHRT_MAX)
                  {
                     /* error - overflow */
                     cast_max = SHRT_MAX;
                  }
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     cast_min, cast_max);
               }
               break;
               default:
                  /* error invalid type */
                  dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                           Generic, out_names[rr]);
               break;
            }
         } 
         break;
         case HDET_DET_X:   
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0,
                             DS_HRC_DET_TUNITS, "Focal plane X");
            if ((inp_p->hrc_system == HRC_IMG_SYS) || 
                (inp_p->hrc_system == HSI_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                        DS_HRCI_DETX_TLMIN, DS_HRCI_DETX_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                        (float) DS_HRCI_DETX_TLMIN, 
                        (float) DS_HRCI_DETX_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long) DS_HRCI_DETX_TLMIN, 
                        (long) DS_HRCI_DETX_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short) DS_HRCI_DETX_TLMIN,
                           cast_max = (short) DS_HRCI_DETX_TLMAX;
 
                     if (DS_HRCI_DETX_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCI_DETX_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else if ((inp_p->hrc_system == HRC_SPC_SYS) || 
                     (inp_p->hrc_system == HRC_SPC_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                        DS_HRCS_DETX_TLMIN, DS_HRCS_DETX_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                        (float) DS_HRCS_DETX_TLMIN, 
                        (float) DS_HRCS_DETX_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long) DS_HRCS_DETX_TLMIN, 
                        (long) DS_HRCS_DETX_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short)DS_HRCS_DETX_TLMIN,
                           cast_max = (short)DS_HRCS_DETX_TLMAX;
 
                     if (DS_HRCS_DETX_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCS_DETX_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else
            {
               /* error - unknown system type */
               dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                        Generic, out_names[rr]);
            }
         break;
         case HDET_DET_Y:   
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0,
                             DS_HRC_DET_TUNITS, "Focal plane Y");
            if ((inp_p->hrc_system == HRC_IMG_SYS) || 
                (inp_p->hrc_system == HSI_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                        DS_HRCI_DETY_TLMIN, DS_HRCI_DETY_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                        (float) DS_HRCI_DETY_TLMIN, 
                        (float) DS_HRCI_DETY_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long) DS_HRCI_DETY_TLMIN, (long) DS_HRCI_DETY_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short) DS_HRCI_DETY_TLMIN,
                           cast_max = (short) DS_HRCI_DETY_TLMAX;

                     if (DS_HRCI_DETY_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCI_DETY_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else if ((inp_p->hrc_system == HRC_SPC_SYS) || 
                     (inp_p->hrc_system == HRC_SPC_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                       DS_HRCS_DETY_TLMIN,DS_HRCS_DETY_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                       (float) DS_HRCS_DETY_TLMIN,
                       (float) DS_HRCS_DETY_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long)DS_HRCS_DETY_TLMIN, (long)DS_HRCS_DETY_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short)DS_HRCS_DETY_TLMIN,
                           cast_max = (short)DS_HRCS_DETY_TLMAX;

                     if (DS_HRCS_DETY_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCS_DETY_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                                            cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else
            {
               /* error - unknown system type */
               dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                        Generic, out_names[rr]);
            }
         break;
         case HDET_SKY_X:   
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0,
                             DS_HRC_SKY_TUNITS, "Sky X");
            if ((inp_p->hrc_system == HRC_IMG_SYS) || 
                (inp_p->hrc_system == HSI_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                        DS_HRCI_X_TLMIN, DS_HRCI_X_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                        (float) DS_HRCI_X_TLMIN, (float) DS_HRCI_X_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long) DS_HRCI_X_TLMIN, (long) DS_HRCI_X_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short) DS_HRCI_X_TLMIN,
                           cast_max = (short) DS_HRCI_X_TLMAX;
 
                     if (DS_HRCI_X_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCI_X_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else if ((inp_p->hrc_system == HRC_SPC_SYS) || 
                     (inp_p->hrc_system == HRC_SPC_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                       DS_HRCS_X_TLMIN,DS_HRCS_X_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                        (float) DS_HRCS_X_TLMIN, (float) DS_HRCS_X_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long)DS_HRCS_X_TLMIN, (long)DS_HRCS_X_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short)DS_HRCS_X_TLMIN,
                           cast_max = (short)DS_HRCS_X_TLMAX;
 
                     if (DS_HRCS_X_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCS_X_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else
            {
               /* error - unknown system type */
               dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                        Generic, out_names[rr]);
            }
         break;
         case HDET_SKY_Y:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, 
                             DS_HRC_SKY_TUNITS, "Sky Y");
            if ((inp_p->hrc_system == HRC_IMG_SYS) || 
                (inp_p->hrc_system == HSI_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                        DS_HRCI_Y_TLMIN, DS_HRCI_Y_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                        (float) DS_HRCI_Y_TLMIN, (float) DS_HRCI_Y_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long) DS_HRCI_Y_TLMIN, (long) DS_HRCI_Y_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short) DS_HRCI_Y_TLMIN,
                           cast_max = (short) DS_HRCI_Y_TLMAX;
    
                     if (DS_HRCI_Y_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCI_Y_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else if ((inp_p->hrc_system == HRC_SPC_SYS) || 
                     (inp_p->hrc_system == HRC_SPC_IMG_SYS))
            {
               switch(evtout_p->types[rr])
               {
                  case dmDOUBLE:
                     dmDescriptorSetRange_d(evtout_p->desc[rr],
                       DS_HRCS_Y_TLMIN,DS_HRCS_Y_TLMAX);
                  break;
                  case dmFLOAT:
                     dmDescriptorSetRange_f(evtout_p->desc[rr],
                       (float) DS_HRCS_Y_TLMIN, (float) DS_HRCS_Y_TLMAX);
                  break;
                  case dmLONG:
                     dmDescriptorSetRange_l(evtout_p->desc[rr],
                        (long)DS_HRCS_Y_TLMIN, (long)DS_HRCS_Y_TLMAX);
                  break;
                  case dmSHORT: 
                  {
                     short cast_min = (short)DS_HRCS_Y_TLMIN,
                           cast_max = (short)DS_HRCS_Y_TLMAX;
 
                     if (DS_HRCS_Y_TLMIN < SHRT_MIN)
                     {
                        /* error - underflow */
                        cast_min = SHRT_MIN;
                     }
                     if (DS_HRCS_Y_TLMAX > SHRT_MAX)
                     {
                        /* error - overflow */
                        cast_max = SHRT_MAX;
                     }
                     dmDescriptorSetRange_s(evtout_p->desc[rr],
                        cast_min, cast_max);
                  }
                  break;
                  default:
                     /* error invalid type */
                     dsErrAdd(hpe_err_p, dsAPEDMDATATYPEERR, Individual, 
                              Generic, out_names[rr]);
                  break;
               }
            }
            else
            {
               /* error - unknown system type */
               dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual, 
                        Generic, out_names[rr]);
            }
         break; 
         case HDET_CP_X:         /*crsu*/
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                              evtout_p->types[rr],  0, "", 
                              "Coarse position U axis");
         break;
         case HDET_CP_Y:         /*crsv*/
            evtout_p->desc[rr] =
               dmColumnCreate(evtout_p->extension, out_names[rr],
                              evtout_p->types[rr],  0, "", 
                              "Coarse position V axis");
         break;
         case HDET_AX_1:         /*au1*/
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", 
                             "U axis ADC 1");
         break; 
         case HDET_AX_2:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", 
                             "U axis ADC 2");
         break;
         case HDET_AX_3:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", "U axis ADC 3");
         break;
         case HDET_AY_1:         /* av1 */
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", "V axis ADC 1");
         break; 
         case HDET_AY_2:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", "V axis ADC 2");
         break; 
         case HDET_AY_3:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", "V axis ADC 3");
         break; 
         case HDET_SUMAMPS:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", 
                             "Sum of all amp readouts");
            dmDescriptorSetRange_s(evtout_p->desc[rr],
               DS_HRC_SUMAMPS_TLMIN, DS_HRC_SUMAMPS_TLMAX);
         break; 
         case HDET_CHIP_ID:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", "Chip ID");
            switch (inp_p->hrc_system)
            {
               case HRC_IMG_SYS:  /* fallthrough intended */ 
               case HSI_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCI_CHIP_ID_TLMIN, DS_HRCI_CHIP_ID_TLMAX);
               break;
               case HRC_SPC_SYS:  /* fallthrough intended */
               case HRC_SPC_IMG_SYS:
                  dmDescriptorSetRange_s(evtout_p->desc[rr],
                     DS_HRCS_CHIP_ID_TLMIN, DS_HRCS_CHIP_ID_TLMAX);
               break;
               default:
                  /* error - unknown system type */
                  dsErrAdd(hpe_err_p, dsAPEUNKNOWNSYSERR, Individual,
                           Generic, out_names[rr]);
               break;
            }
         break; 
         case HDET_AMP_SF:
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", 
                             "Amplitude scale factor");
         break; 
         case HDET_STATUS:
            if (evtout_p->types[rr] == dmBIT)
            { 
               evtout_p->desc[rr] =
                  dmColumnCreateArray(evtout_p->extension, out_names[rr],
                                      evtout_p->types[rr],  0, "", 
                                      "Event status bits", 4);
            }
            else
            {
               evtout_p->desc[rr] =
                  dmColumnCreate(evtout_p->extension, out_names[rr],
                                 evtout_p->types[rr],  0, "", 
                                 "Event status bits");
            } 
         break; 
         default :
            evtout_p->desc[rr] =
              dmColumnCreate(evtout_p->extension, out_names[rr],
                             evtout_p->types[rr],  0, "", "");
         break; 
      }
   }
} /*end: hrc_setup_columns */
