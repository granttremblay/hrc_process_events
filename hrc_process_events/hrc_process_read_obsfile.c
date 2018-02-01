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

/*H**************************************************************************
 
* FILE NAME: hrc_process_read_obsfile.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: 
 
* NOTES:

*JCC(8/15/00) - display an error when paramopen returns NULL.
*JCC(8/2002)-implement amp_sf correction :
*JCC(11/2002)- get inp_p->ra_nom/dec_nom from obsfile (inp_p->obs_info_p).

(Algorithm of getting RA_NOM, DEC_NOM 
  -set inp_p->ra_nom=NOM_NOT_FOUND=-9999.0; 
  -if (obsfile != NONE)
      - get NOM from obsfile  ;
      - if NOM found in obsfile ;
            inp_p->ra_nom=obsfile's ;
   else
      - do nothing ;
  -If [ find NOM in evt infile  ( =evt_ra_nom ) ]
      - save the value as  'evt_ra_nom' ;
      - if ( inp_p->ra_nom== NOM_NOT_FOUND) -> not found in obsfile or obsfile==NONE
          inp_p->ra_nom = evt_ra_nom ;
      - else
          if (evt_ra_nom != inp_p->ra_nom )
             WARNING:  obsfile and infile have different NOM values and
                       will take it from obsfile 
      - endif
   Else
      - if (inp_p->ra_nom==NOM_NOT_FOUND) -> obsfile!=NONE, missing in obsfile+infile;
                                          -> obsfile==NONE, missing in infile;
         ERROR: Require RA/DEC_NOM keywords in obsfile or evtfile.
)
*JCC(3/2003) - replace old hdrlib with new one  ( new getHdr )
*              call hdrGetKeyValue ; hdrFindKey ;
*
*JCC(4/2003)-new tapRing spec ;
*H**************************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 


void hrc_process_read_obsfile(
   INPUT_PARMS_P_T  inp_p,     /* I/O - input parameters/data        */
   dsErrList*       err_p)     /* O   - error list pointer           */
{

   char  *tmp ;

   /* check if obs.par provided */ 
   /* if  ( obsfile != "none" )   */
   if (ds_strcmp_cis(inp_p->obsfile, "none") != 0)
   {
   /* ---------------------------------------------------------------------
    * JCC(8/00)-display error for a compressed file.
    * ---------------------------------------------------------------------*/
    /* if obs file exists, open it */
    paramfile  *obs_par_file;
    obs_par_file = paramopen(inp_p->obsfile, NULL, 0, "rH");

    if (obs_par_file != NULL )
    {
      /* don't clobber event or primary header */
      if (inp_p->obs_info_p == NULL) 
      {
         /* (3/2003) - no need to allocate memory for obs_info_p for new hdrlib */
         /* (3/2003) - getHdr will allocate memory for obs_info_p */

            if (!ds_check_access(inp_p->obsfile,"R"))
            {
               /* ---------------------------------------------
                * (8/2002)- amp_sf correction :
                *
                * RANGE_SWITCH_LEVEL  will be got either from 
                * obspar file (for AP run)   or from event lev1 
                * ( for user's run ),  otherwise, give ERROR.
                *   
                * obspar supersedes the event file.
                * ---------------------------------------------*/
               if (inp_p->do_amp_sf_cor == TRUE)
               {
                  /*-----------------------------------------------------
                   * pgets : for data type 'short' ; 
                   *
                   * if not found, the ERROR message from pgets() ;
                   *  obs0a.par:  parameter not found : range_switch_level ;
                   *----------------------------------------------------*/
                   if (inp_p->get_range_switch_level==FALSE) /*not get yet*/
                   {
                      /* output keyword 'RANGELEV' */
                      if(paccess(obs_par_file, "range_switch_level"))
                      {
                          inp_p->range_switch_level = pgets(obs_par_file, 
                                                  "range_switch_level");
                          inp_p->get_range_switch_level = TRUE ;
                          inp_p->AMPSFCOR = TRUE ;      /* for outfile */

                          /* always get   obs_r_s_l */
                          inp_p->obs_r_s_l = inp_p->range_switch_level ;
                      }
                      else
                      {
                         dsErrAdd(err_p, dsGENERICERR, Individual,Custom,
                            " WARNING:  range_switch_level  not found in %s.", 
                            inp_p->obsfile );
                         inp_p->get_range_switch_level = FALSE ;
                         inp_p->AMPSFCOR = FALSE ;      /* for outfile */
                      }
                   }
               }    
               else
               {
                   inp_p->get_range_switch_level = FALSE ;
                   inp_p->range_switch_level = 0 ;   /* default; */ 
                   inp_p->AMPSFCOR = FALSE ;       

                   /* always get    obs_r_s_l */
                   if(paccess(obs_par_file, "range_switch_level"))
                   {
                      inp_p->obs_r_s_l = pgets(obs_par_file,
                                               "range_switch_level");
                   }
               }
               /* ------ end : (8/2002) ----------------*/

              /*-----------------------------------------------------
               * (4/2003) get width_threshold from obs.par
               *----------------------------------------------------*/
               if (paccess(obs_par_file, "width_threshold"))
               {
                   inp_p->obs_widthres = pgets(obs_par_file, "width_threshold");
               }
               else
               {
                   /* dsErrAdd(err_p, dsGENERICERR, Individual,Custom, " WARNING:  width_threshold not found in %s.", inp_p->obsfile ); */
               }


              /*---------------------------
               * load event header values 
               *---------------------------*/
               /* (2/2003) - replace old get_event_hdr with new getHdr */
               inp_p->obs_info_p = getHdr(obs_par_file, hdrPAR_FILE );

               /* (2/2003)- get key values for TSTART TSTOP SIM_X/Y/Z */
               hdrGetKeyValue_d( inp_p->obs_info_p, "TSTART", &(inp_p->obs_tstart)) ;
               hdrGetKeyValue_d( inp_p->obs_info_p, "TSTOP", &(inp_p->obs_tstop)) ;
               hdrGetKeyValue_d( inp_p->obs_info_p, "SIM_X", &(inp_p->sim_x)) ;
               hdrGetKeyValue_d( inp_p->obs_info_p, "SIM_Y", &(inp_p->sim_y)) ;
               hdrGetKeyValue_d( inp_p->obs_info_p, "SIM_Z", &(inp_p->sim_z)) ;

               /* (2/2003)- get key values for RA_NOM/DEC_NOM if they are found */
               if ( hdrFindKey(inp_p->obs_info_p, "RA_NOM") == hdrFOUND )
                  hdrGetKeyValue_d( inp_p->obs_info_p, "RA_NOM", &(inp_p->ra_nom)) ;
               if ( hdrFindKey(inp_p->obs_info_p, "DEC_NOM") == hdrFOUND )
                  hdrGetKeyValue_d( inp_p->obs_info_p, "DEC_NOM", &(inp_p->dec_nom)) ;

               inp_p->use_obs = 1; 

               /* (2/2003)- get values of DATAMODE DETNAM TELESCOP */ 
               /*  DS_SZ_KEYWORD = HDR_SZ = 511 */
               memset( inp_p->datamode, 0, DS_SZ_KEYWORD);
               memset( inp_p->instrume, 0, DS_SZ_KEYWORD);
               memset( inp_p->telescop, 0, DS_SZ_KEYWORD);
               hdrGetKeyValue_c( inp_p->obs_info_p, "DATAMODE", &tmp );
               strcpy( inp_p->datamode, tmp ) ;
               inp_p->next_in_line =
                  ((ds_strcmp_cis(inp_p->datamode, NEXT_IN_LINE_VAL) == 0) ||
                   (ds_strcmp_cis(inp_p->datamode, NEXT_IN_LINE_VAL2) == 0));

               hdrGetKeyValue_c( inp_p->obs_info_p, "DETNAM", &tmp ) ;
               strcpy( inp_p->instrume, tmp ) ;

               hdrGetKeyValue_c( inp_p->obs_info_p, "TELESCOP", &tmp ) ;
               strcpy( inp_p->telescop, tmp ) ;
               if ((ds_strcmp_cis(inp_p->telescop, TELESCOP_XRCF_KEY) == 0) ||
                   (ds_strcmp_cis(inp_p->telescop, TELESCOP_XRCF_KEY2) == 0))
               {
                  inp_p->processing = HRC_PROC_XRCF;
               }


               /* load primary header info with get_prim_hdr (removed, 2/2003)*/ 
               /* (2/2003) - no need to call getHdr twice ; */

               paramclose(obs_par_file);
            }  /* end :  if (!ds_check_access(inp_p->obsfile,"R")) */
      }  /* end: if(inp_p->obs_info_p==NULL) */
      else 
      {
         /* error- obs.par already read */ 
      }  /* end: if (inp_p->obs_info_p==NULL) */
    }  /* end: if (obs_par_file != NULL ) */
    else
    {
        dsErrAdd(err_p, dsREADFILEFERR, Individual, Generic, inp_p->obsfile);
    }  /* end: if (obs_par_file != NULL ) */
   } /*end:   if (ds_strcmp_cis(inp_p->obsfile,"none")!=0) */
}
