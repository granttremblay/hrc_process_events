/*                                                                
**  Copyright (C) 1996-2009  Smithsonian Astrophysical Observatory 
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
 
* FILE NAME: hrc_process_events_functions.c
*
* DEVELOPEMENT: tools
*
* DESCRIPTION: This routine contains several functions used by 
* hrc_process_events. The modules in this file are:
*
*    load_input_parameters()
*    set_up_mirror() 
*    convert_coord_type_to_char()
*    adjust_output_eventdef()
*
* For details on the functionality of any of the above listed modules, 
* please see the 'description' comment preceding the specific module's
* source code.
*
*
* REVISION HISTORY:
* wmclaugh(8/19/96) - First Version.
* JCC(5/1/00)-update load_input_parameters to read 'tapfile' (ie. tap ring 
*             coefficients filename )
* JCC(2/2002)- pass geompar to the pixlib call (load_input_parameters)
* JCC(8/2002)- add do_amp_sf_cor, ampsfcorfile, range_switch_level,
*              get_range_switch_level and AMPSFCOR to load_input_parameters;
*            - initialize AMPSFCOR,range_switch_level 
*              in load_input_parameters ;
*            - inp_p->range_switch_level is 'range_switch_level' from obspar
*              or RANGELEV from evt1 header file.
* (1/2009)- Add gdropfile parameter for hrcS 3dim gain:  ( obsolete 10/2009 )
*10/2009 - remove gdropfile from hpe.par
*H***********************************************************************/

#include <float.h> 

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

#include "parameter.h"

/*************************************************************************
 
* DESCRIPTION:

  The routine set_up_mirror() is called by hrc_process_events to set the 
  staging/fam data needed by the pixel library to perform coordinate 
  transformations. The routine attempts to load keyword values provided in 
  the input file headers. If the keywords do not exist in the file headers,
  default values of 0 are used. 
 
  Sim keywords are obtained as follows. First, if an obs.par file was
  provided, those sim_x/y/z keyword values are used. If an obs.par file
  is not utilized, the routine attempts to read the keywords from the
  principal extension of the first input event file. If the keywords do
  not exist in the principal extension and no obs.par file was provided,
  the default pixlib configuration is utilized. 
  
*************************************************************************/

void set_up_mirror(
   dmBlock*          principal, /* I - pntr to input event file handle  */
   INPUT_PARMS_P_T   inp_p,     /* I - input parameter structure        */ 
   char*             evtfile,   /* I - name of input event file         */
   ALIGNMENT_REC_P_T aln_p,     /* O - pointer to alignment record      */
   dsErrList* err_p)            /* O - pointer to error message stack   */
{
   if ((dmKeyRead_d(principal, STG_X_KEY, &aln_p->pos[STG_AIMPOINT_X]) == 
       NULL) ||
       (dmKeyRead_d(principal, STG_Y_KEY, &aln_p->pos[STG_AIMPOINT_Y]) == 
       NULL) ||
       (dmKeyRead_d(principal, STG_Z_KEY, &aln_p->pos[STG_AIMPOINT_Z]) == 
       NULL))
   {
      aln_p->pos[STG_AIMPOINT_X] = aln_p->pos[STG_AIMPOINT_Y] = 0.0;
      aln_p->pos[STG_AIMPOINT_Z] = 0.0;
      if (inp_p->processing != HRC_PROC_FLIGHT)
      {
         dsErrAdd(err_p, dsHPESTAGEDNEERR, Individual, Generic, evtfile);
      }
   }

   if ((dmKeyRead_d(principal, STG_ANG1_KEY, &aln_p->theta[STG_ANG1_NDX]) ==
       NULL) ||
       (dmKeyRead_d(principal, STG_ANG2_KEY, &aln_p->theta[STG_ANG2_NDX]) ==
       NULL) ||
       (dmKeyRead_d(principal, STG_ANG3_KEY, &aln_p->theta[STG_ANG3_NDX]) ==
       NULL))
   {
      aln_p->theta[STG_ANG1_NDX] = aln_p->theta[STG_ANG2_NDX] = 0.0;
      aln_p->theta[STG_ANG3_NDX] = 0.0;
      if (inp_p->processing != HRC_PROC_FLIGHT)
      {
         dsErrAdd(err_p, dsHPESTAGEANGDNEERR, Individual, Generic, evtfile);
      }
   }

   if ((dmKeyRead_d(principal, HPY_YAW_KEY, &aln_p->hpy[HPY_YAW]) == NULL) ||
       (dmKeyRead_d(principal, HPY_PIT_KEY, &aln_p->hpy[HPY_PIT]) == NULL))
   {
      aln_p->hpy[HPY_YAW] = aln_p->hpy[HPY_PIT] = 0.0;
      if (inp_p->processing != HRC_PROC_FLIGHT)
      {
         dsErrAdd(err_p, dsHPEHPYDNEERR, Individual, Generic, evtfile);
      }
   }
 
   pix_set_detector(inp_p->instrume);

   if (inp_p->processing == HRC_PROC_XRCF)
   {
      /* xrcf- only need sim_z keyword to overide default aimpoint */
      if (inp_p->use_obs)
      {
         /* if an obs.par file was provided, use the sim_z value */ 
         aln_p->lass[SIM_Z_NDX] = inp_p->sim_z; 
	 pix_set_aimpoint_by_value(aln_p->lass);
      }
      else if (dmKeyRead_d(principal, SIM_Z_KEY, &aln_p->lass[SIM_Z_NDX]) 
         != NULL)
      {
         /* if no obs.par, try to get the EVENTS extension keyword */
         pix_set_aimpoint_by_value(aln_p->lass);
      }
      else
      {
         /* get the pixlib default positions */ 
	 pix_set_aimpoint_default();
	 pix_get_aimpoint_by_value(aln_p->lass);

         dsErrAdd(err_p, dsHPESETAIMPOINTERR, Individual, Generic, evtfile);
      }
   }
   else
   {
      /* flight- need all 3 sim keywords to overide default aimpoint */
      if (inp_p->use_obs)
      {
         /* if an obs.par file was provided, use the sim values */
         aln_p->lass[SIM_X_NDX] = inp_p->sim_x;
         aln_p->lass[SIM_Y_NDX] = inp_p->sim_y;
         aln_p->lass[SIM_Z_NDX] = inp_p->sim_z;
         pix_set_aimpoint_by_value(aln_p->lass);
      }
      else if ((dmKeyRead_d(principal, SIM_X_KEY, &aln_p->lass[SIM_X_NDX]) 
          != NULL) &&
          (dmKeyRead_d(principal, SIM_Y_KEY, &aln_p->lass[SIM_Y_NDX]) != 
          NULL) &&
          (dmKeyRead_d(principal, SIM_Z_KEY, &aln_p->lass[SIM_Z_NDX]) != 
          NULL))
      {
         /* if no obs.par, try to get the EVENTS extension keywords */
         pix_set_aimpoint_by_value(aln_p->lass);
      }
      else
      {
         /* get the pixlib default positions */ 
	 pix_set_aimpoint_default();
	 pix_get_aimpoint_by_value(aln_p->lass);

         dsErrAdd(err_p, dsHPESETAIMPOINTERR, Individual, Generic, evtfile);
      }
   }
 
   pix_set_mirror (aln_p->hpy, aln_p->pos, aln_p->theta);

} 


/****************************************************************************
*
* DESCRIPTION:
*
* The routine load_input_parameters() is called upon by hrc_process_events to 
* obtain input parameters from the hrc_process_events.par parameter file. The 
* routine does not validate the input. It merely copies the values from the 
* parameter file into the appropriate field of a data structure of type 
* INPUT_PARMS_T. If the parameter file does not contain an expected param,
* an error message is added to the error list.
*
* JCC(2/2002) - pass geompar to the pixlib call
* JCC(8/2002) - get do_amp_sf_cor, ampsfcorfile ;
*               initialize range_switch_levelt, get_range_switch_level ; 
*             - initialize use_obs (1=access_obsfile_successfully) ;
*JCC(4/2003)-update load_input_parameters to intialize inp_p; 
*       - rearrang the order of calling clgstr;
****************************************************************************/

void load_input_parameters(
   INPUT_PARMS_P_T inp_p,     /* O - input parameters for the program      */
   dsErrList*      err_p)     /* O - pntr to error message stack           */
{
   /* initialize fields */
   memset(inp_p, 0, sizeof(INPUT_PARMS_T)); 
   inp_p->processing = HRC_PROC_FLIGHT;
   inp_p->evt_tstart = DBL_MAX;
   inp_p->evt_tstop = DBL_MIN;
   inp_p->gain_cdelt[0] = GAIN_DEFAULT_CDELT_X; 
   inp_p->gain_cdelt[1] = GAIN_DEFAULT_CDELT_Y; 

   /*-----------------------------------------------------------
    * JCC(8/2002) initialize AMPSFCOR,range_switch_level(=RANGELEV) 
    *              for amp_sf_corrections 
    *----------------------------------------------------------*/
   inp_p->use_obs = 0 ; 
   inp_p->get_range_switch_level = FALSE ;
   inp_p->range_switch_level = 0 ;
   inp_p->AMPSFCOR = FALSE ;  /* AMPSFCOR from obsfile or evt1 infile*/
   inp_p->evt_AMPSFCOR = FALSE;  /* AMPSFCOR from evt1 infile */
   inp_p->match_range_switch_level=FALSE; /*find the matched r_s_l in caldb? */ 

   inp_p->obs_r_s_l = NEG_9999 ;  /* range_switch_level should be >=0.0 */
   inp_p->evt_r_s_l = NEG_9999 ;  /* RANGELEV from evt infile */

   inp_p->obs_widthres = NEG_9999 ;  /* width_thresold from obs.par */
   inp_p->evt_widthres = NEG_9999 ;  /* WIDTHRES from evt1 file */

   inp_p->ra_nom  = NOM_NOT_FOUND ;  /* -9999.0 */
   inp_p->dec_nom = NOM_NOT_FOUND ;  /* -9999.0 */

   /*-----------------------------------------------------------*/

   /* check if parameter exists before attempting to open it */
   if (paccess(PFFile, "infile"))
   {
      clgstr("infile", inp_p->stack_in, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "infile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "outfile"))
   {
      clgstr("outfile", inp_p->outfile, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "outfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "badpixfile"))
   {
      clgstr("badpixfile", inp_p->badpixfile, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "badpixfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "acaofffile"))
   {
      clgstr("acaofffile", inp_p->asp_file, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "acaofffile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "alignmentfile"))
   {
      clgstr("alignmentfile", inp_p->align_file, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "alignmentfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "obsfile"))
   {
      clgstr("obsfile", inp_p->obsfile, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "obsfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "geompar"))
   {
      clgstr("geompar", inp_p->geompar, DS_SZ_PATHNAME);
   }
   else
   {
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "geompar", "hrc_process_events.par");
   }
   if (paccess(PFFile, "do_ratio"))
   {
      inp_p->do_ratio  = clgetb("do_ratio");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "do_ratio", "hrc_process_events.par");
   }
   if (paccess(PFFile, "do_amp_sf_cor"))
   {
      inp_p->do_amp_sf_cor = clgetb("do_amp_sf_cor");
      inp_p->get_range_switch_level=FALSE; /*we've NOT got RANGE_SWITCH_LEVEL from either obs.par or event1 file */
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "do_amp_sf_cor", "hrc_process_events.par");
   }
   if (paccess(PFFile, "gainfile"))
   {
      clgstr("gainfile", inp_p->gain_file, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "gainfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "ADCfile"))
   {
      clgstr("ADCfile", inp_p->adc_file, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "ADCfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "degapfile"))
   {
      clgstr("degapfile", inp_p->degap_file, HRC_DEGAP_FILE_LEN);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "degapfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "hypfile"))
   {
      clgstr("hypfile", inp_p->hypfile, HRC_DEGAP_FILE_LEN);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "hypfile", "hrc_process_events.par");
   }
   if (inp_p->do_amp_sf_cor)
   {
      if (paccess(PFFile, "ampsfcorfile"))
      {
         clgstr("ampsfcorfile", inp_p->ampsfcorfile, HRC_DEGAP_FILE_LEN);
      }
      else
      {
         /* parameter does not exist- push error message onto stack */
         dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
                  "ampsfcorfile", "hrc_process_events.par");
      }
   }
   if (paccess(PFFile, "tapfile"))
   {
      clgstr("tapfile", inp_p->tapfile, HRC_DEGAP_FILE_LEN);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "tapfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "ampsatfile"))
   {
      clgstr("ampsatfile", inp_p->ampsatfile, HRC_DEGAP_FILE_LEN);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "ampsatfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "evtflatfile"))
   {
      clgstr("evtflatfile", inp_p->ampflatfile, HRC_DEGAP_FILE_LEN);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "evtflatfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "badfile"))
   {
      clgstr("badfile", inp_p->badfile, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "badfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "logfile"))
   {
      clgstr("logfile", inp_p->logfile, DS_SZ_PATHNAME);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "logfile", "hrc_process_events.par");
   }
   if (paccess(PFFile, "instrume"))
   {
      clgstr("instrume", inp_p->instrume, DS_SZ_KEYWORD);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "instrume", "hrc_process_events.par");
   }
   if (paccess(PFFile, "eventdef"))
   {
      clgstr("eventdef", inp_p->outcols, DS_SZ_COMMAND);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "eventdef", "hrc_process_events.par");
   }
   if (paccess(PFFile, "badeventdef"))
   {
      clgstr("badeventdef", inp_p->badoutcols, DS_SZ_COMMAND);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "badeventdef", "hrc_process_events.par");
   }
   if (paccess(PFFile, "grid_ratio"))
   {
      inp_p->grid_ratio = clgetd("grid_ratio");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "grid_ratio", "hrc_process_events.par");
   }
   if (paccess(PFFile, "pha_ratio"))
   {
      inp_p->pha_ratio = clgetd("pha_ratio");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "pha_ratio", "hrc_process_events.par");
   }
   if (paccess(PFFile, "wire_charge"))
   {
      inp_p->wire_charge = clgeti("wire_charge");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "wire_charge", "hrc_process_events.par");
   }
   if (paccess(PFFile, "cfu1"))
   {
      inp_p->cf[HDET_PLANE_X][HDET_1ST_ORD_CF] = clgetd("cfu1");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "cfu1", "hrc_process_events.par");
   }
   if (paccess(PFFile, "cfu2"))
   {
      inp_p->cf[HDET_PLANE_X][HDET_2ND_ORD_CF] = clgetd("cfu2");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "cfu2", "hrc_process_events.par");
   }
   if (paccess(PFFile, "cfv1"))
   {
      inp_p->cf[HDET_PLANE_Y][HDET_1ST_ORD_CF] = clgetd("cfv1");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "cfv1", "hrc_process_events.par");
   }
   if (paccess(PFFile, "cfv2"))
   {
      inp_p->cf[HDET_PLANE_Y][HDET_2ND_ORD_CF] = clgetd("cfv2");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "cfv2", "hrc_process_events.par");
   }
   if (paccess(PFFile, "time_offset"))
   {
      inp_p->time_offset = clgetd("time_offset");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "time_offset", "hrc_process_events.par");
   }
   if (paccess(PFFile, "amp_gain"))
   {
      inp_p->amp_gain = clgetd("amp_gain");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "amp_gain", "hrc_process_events.par");
   }
   if (paccess(PFFile, "rand_seed"))
   {
      inp_p->rand_seed = clgeti("rand_seed");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "rand_seed", "hrc_process_events.par");
   }
   if (paccess(PFFile, "rand_pix_size"))
   {
      inp_p->randpixsize = clgetr("rand_pix_size");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "rand_pix_size", "hrc_process_events.par");
   }
   if (paccess(PFFile, "tstart"))
   {
      clgstr("tstart", inp_p->time_start, DS_SZ_KEYWORD);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "tstart", "hrc_process_events.par");
   }
   if (paccess(PFFile, "tstop"))
   {
      clgstr("tstop", inp_p->time_stop, DS_SZ_KEYWORD);
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "tstop", "hrc_process_events.par");
   }

   if (paccess(PFFile, "start"))
   {
      clgstr("start", inp_p->start_coord, DS_SZ_KEYWORD);
      string_to_lowercase(inp_p->start_coord, DS_SZ_KEYWORD); 
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "start", "hrc_process_events.par");
   }
   if (paccess(PFFile, "stop"))
   {
      clgstr("stop", inp_p->stop_coord, DS_SZ_KEYWORD);
      string_to_lowercase(inp_p->stop_coord, DS_SZ_KEYWORD); 
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "stop", "hrc_process_events.par");
   }
   if (paccess(PFFile, "clobber"))
   {
      inp_p->clobber  = clgetb("clobber");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "clobber", "hrc_process_events.par");
   }
   if (paccess(PFFile, "verbose"))
   {
      inp_p->debug = clgeti("verbose");
   }
   else
   {
      /* parameter does not exist- push error message onto stack */
      dsErrAdd(err_p, dsFINDPARAMFERR, Individual, Generic,
               "verbose", "hrc_process_events.par");
   }

   inp_p->do_ADC = ((ds_strcmp_cis(inp_p->adc_file, "NONE") != 0) && 
                     (strcmp(inp_p->adc_file, "\0") != 0)); 
 
   inp_p->do_pi = ((ds_strcmp_cis(inp_p->gain_file, "NONE") != 0) && 
                   (strcmp(inp_p->gain_file, "\0") != 0)); 
}




/****************************************************************************
 
* FUNCTION NAME: hrc_process_evt_file_cleanup()
 
* DESCRIPTION:
 
  The routine hrc_process_evt_file_cleanup() is called by hrc_process_events
  deallocate the memory used by the data model and close any open file
  descriptors.
 
  The function does not pass a return value back to the calling routine.
 
****************************************************************************/
 
void hrc_process_evt_file_cleanup(EVENT_SETUP_P_T evt_set_p)
{
   if ((evt_set_p->primary != evt_set_p->extension) &&
       (evt_set_p->primary != NULL))
   {
      /* close the primary header extension  */
      dmBlockClose(evt_set_p->primary);
      evt_set_p->primary = NULL;
   }
   if (evt_set_p->extension != NULL)
   {
      /* close the event file */
      dmDatasetTableClose(evt_set_p->extension);
      evt_set_p->extension = NULL;
      evt_set_p->dataset = NULL;
      evt_set_p->primary = NULL;
   }
   if (evt_set_p->desc != NULL)
   {
      /* free dynamic memory for attributes */
      free(evt_set_p->desc);
      evt_set_p->desc = NULL;
   }
   if (evt_set_p->mapping != NULL)
   {
      /* free dynamic memory for mapping */
      free(evt_set_p->mapping);
      evt_set_p->mapping = NULL;
   }
   if (evt_set_p->types != NULL)
   {
      /* free dynamic memory for types data */
      free(evt_set_p->types);
      evt_set_p->types = NULL;
   }
   if (evt_set_p->dim != NULL)
   {
      /* free dynamic memory for dimension data */
      free(evt_set_p->dim);
      evt_set_p->dim = NULL; 
   }
}
