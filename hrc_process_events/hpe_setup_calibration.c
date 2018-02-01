/*                                                                
**  Copyright (C) 1996-2009,2010,2012  Smithsonian Astrophysical Observatory 
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

/************************************************************************
 *JCC(8/2002)-update access_caldb() to retrieve the header keyword 
 *            RANGELEV (in obs.par,it's 'range_switch_level') from event 
 *            lev1 input file.
 *JCC(4/2003)-new tapRing spec; 
 *   - get WIDTHRES from evt1 (in obs par, it's "width_threshold");
 *   - width_threshold :  obs par overwrites evt1 ; give warning 
 *     if not found in evt1 and obs.par
 *(12/2004)- let obspar's detnam override evtfile's (calOverride).
 *(5/2005)-add GAPLOOKUP to get the new #3 degap file
 *(8/2005)-add calQUERY_EXPRESSION for #3 degap file to avoid warning of
 *        (4||2 CALDB files found. Using the first )
 *
 *(5/2007)-new caldb4.
 *   If use_obs=1, get hdr info from obs.par ;
 *   Else   get hdr info from 1st input evt file.
 *(9/2008)- remove calSetParam, calSetMatchMode.
 *(10/2008)- replace ds_strncmp_cis w/ ds_strcmp_cis for caldb4.
 *(1/2009)-add caldb product 'GDROP' for gdropfile. (obsolete 10/2009)
10/2009- dph new hrcS gain table doesn't need gdropfile
         and the caldb4 codename is T_GMAP.
4/2010-remove old caldb.
4/2012-bug 13198: gmap not-found warning.
 ************************************************************************/
#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

void access_caldb(
		  INPUT_PARMS_P_T inp_p,
		  dsErrList *all_err_p
		  )
{
  dmBlock *inBlock;

  Stack evtstack;
  char *evtfile = NULL;   /* 1st file on stack */
  
  dmDescriptor *dmkey=NULL, *dmkey2=NULL, *dmkey3=NULL;
  dmDescriptor *dmkey4=NULL ;

  
  evtstack = stk_build(inp_p->stack_in);
  evtfile  = stk_read_num( evtstack, 1 ) ;

  if ( stk_count(evtstack) >= 1 )
    {
      inBlock = dmTableOpen( evtfile );
      
      if (inBlock != NULL )
	{

   if (inp_p->use_obs==1)
      inp_p->caldb4_hdr = inp_p->obs_info_p;           /*getHdr from obs.par*/
   else
      inp_p->caldb4_hdr = getHdr(inBlock, hdrDM_FILE); /*getHdr from 1st evtin block*/

          /*-------------------------------------------------------
           * (8/2002)-retrieve the header keyword RANGELEV from
           *          the 1st evt1 infile.
           *   (RANGELEV in event header = range_switch_level in obspar) 
           *------------------------------------------------------*/
          if  (( inp_p->do_amp_sf_cor == TRUE ) &&
               ( inp_p->get_range_switch_level == FALSE ) )
          {
              dmkey = dmKeyOpen(inBlock, "RANGELEV");
              if ( dmkey != NULL )
              {
                 inp_p->range_switch_level = dmGetScalar_s( dmkey ) ;
                 inp_p->get_range_switch_level = TRUE ; 

                 inp_p->AMPSFCOR = TRUE ;    /* will be in output header */
              }
              else
              {
                 inp_p->AMPSFCOR = FALSE ;    /* will be in output header */
                 inp_p->get_range_switch_level = FALSE ; 

                 /* give ERROR */ 
                 dsErrAdd(all_err_p, dsGENERICERR, Individual,Custom,
                 " WARNING :  RANGELEV keyword is missing in %s.", evtfile);
              }
          }
          /*-------------------------------------------------------
           * (8/2002)-get evt_AMPSFCOR from evt1 header; default=FALSE;
           *-------------------------------------------------------*/
           dmkey2 = dmKeyOpen(inBlock, "AMPSFCOR") ; 
           if ( dmkey2 != NULL )
                inp_p->evt_AMPSFCOR = dmGetScalar_q ( dmkey2 ) ;
          /*-------------------------------------------------------
           * (8/2002)-get evt_r_s_l from evt1 header; default=-9999;
           * (do_amp_sf_cor=no, write out r_s_l that comes from obsfile 
           *  or from evt infile; obsfile supersedes evt;
           *-------------------------------------------------------*/
           dmkey3 = dmKeyOpen(inBlock, "RANGELEV") ;
           /*---- missing r_s_l in obsfile ----*/
           if (( dmkey3 != NULL )&&( inp_p->obs_r_s_l < 0 )) 
                 inp_p->evt_r_s_l = dmGetScalar_s( dmkey3 ) ;
           /*----------------
            * end (8/2002) 
            *----------------*/
          /*------------------------------------------------------------
           *(4/2003)-retrieve the keyword WIDTHRES from 1st evt1 infile.
           * (WIDTHES in event header = width_threshold in obspar)
           *-----------------------------------------------------------*/
           dmkey4 = dmKeyOpen(inBlock, "WIDTHRES") ;
           /*---- width_threshold missing in obsfile ----*/
           if ( inp_p->obs_widthres < 0 )
           {
              if ( dmkey4 != NULL )
              {
                 inp_p->evt_widthres = dmGetScalar_s( dmkey4 ) ;
              }
              else
              {
                 /* dsErrAdd(all_err_p, dsGENERICERR, Individual,Custom, " WARNING : width_threshold/WIDTHRES not found in obs par and evt file.\n"); */
              }
           }
           /* printf("inp_p->obs_widthres=%d\n", inp_p->obs_widthres); */
           /* printf("inp_p->evt_widthres=%d\n", inp_p->evt_widthres); */

          inp_p->hcp = init_caldb_var ( inp_p ) ;

          dsErrCode  Tmp = dsNOERR ;

	  /* ---- Degap ---- */
          Tmp = find_caldb_file( inp_p->degap_file, "GAPLOOKUP", inp_p->hcp ) ;
          if ( Tmp != dsNOERR )
          {
             if ( inp_p->hcp->flg  == INIT_NOT_OK )
                strcpy( inp_p->degap_file, "NONE" ) ;
             else
             {
                dsErrCode  Tmp2 = dsNOERR ;
                Tmp2 = find_caldb_file( inp_p->degap_file, "DEGAP", inp_p->hcp ) ;
                if ( Tmp2 != dsNOERR )
                   strcpy( inp_p->degap_file, "NONE" ) ;
             }
          }
         /* --------------------------------------------------------------- 
          * caldb  gain map  file  :
          * 10/2009 - add T_GMAP for caldb4 file search.
          * ---------------------------------------------------------------*/
          Tmp = find_caldb_file( inp_p->gain_file, "T_GMAP", inp_p->hcp ) ;

         /* 10/2009 - Note :
          * if Tmp=dsNOERR, it could mean : 
          *      [gainfile =(none||" "||"anyfile.fits")]
          *   or [(gainfile="caldb")&&(caldb4 find a file w/ T_GMAP)]
          *
          *   so, don't set the value for gainflag here.
          */
          short gain_warn = 0 ;    /* 4/2012: bug 13198 */
          if ( Tmp != dsNOERR )
          {
             if ( inp_p->hcp->flg  == INIT_NOT_OK )
             {
                gain_warn = 1 ;
                strcpy( inp_p->gain_file, "NONE" ) ;
             }
             else
             {
                dsErrCode  Tmp2 = dsNOERR ;
                Tmp2 = find_caldb_file( inp_p->gain_file, "GMAP", inp_p->hcp);
                if ( Tmp2 != dsNOERR )
                {
                   gain_warn = 1 ;
                   strcpy( inp_p->gain_file, "NONE" ) ;
                }
             }
          }

          if ( gain_warn != 0 )
          {
             err_msg("WARNING: gain map file is not found. set to NONE. \n");
          }

	  /* ---- ADC correction ---- */
          Tmp = find_caldb_file( inp_p->adc_file, "ADC", inp_p->hcp ) ;
          if ( Tmp != dsNOERR )
             strcpy( inp_p->adc_file, "NONE" ) ;

	  /* ---- Hyperbolic test ---- */
          Tmp = find_caldb_file( inp_p->hypfile, "FPTEST", inp_p->hcp ) ;
          if ( Tmp != dsNOERR )
             strcpy( inp_p->hypfile, "NONE" ) ;

	  /* ---- Tap ringing test ---- */
          Tmp = find_caldb_file( inp_p->tapfile, "TAPRINGTEST", inp_p->hcp ) ;
          if ( Tmp != dsNOERR )
             strcpy( inp_p->tapfile, "NONE" ) ;

	  /* ---- Event flatness test ---- */
          Tmp = find_caldb_file( inp_p->ampflatfile, "EFTEST", inp_p->hcp ) ;
          if ( Tmp != dsNOERR )
             strcpy( inp_p->ampflatfile, "NONE" ) ;

	  /* ---- Event saturatio test ---- */
          Tmp = find_caldb_file( inp_p->ampsatfile, "SATTEST", inp_p->hcp ) ;
          if ( Tmp != dsNOERR )
             strcpy( inp_p->ampsatfile, "NONE" ) ;

	  /* ---- AMF_SF Correction ---- */
          Tmp = find_caldb_file( inp_p->ampsfcorfile, "AMP_SF_COR", inp_p->hcp ) ;
          if ( Tmp != dsNOERR )
             strcpy( inp_p->ampsfcorfile, "NONE" ) ;

	  dmTableClose( inBlock );
	} /* end if inblock != NULL */


    } /* end stk count > 1 */

  if (evtstack) stk_close(evtstack);

  return;

} /* end :  access_caldb() */

/*---------------------------------------------------------------------
 * Description : 
 *
 *  find_caldb_file() returns the following Error Status : 
 *
 *  1. dsNOERR :  either myFile is not CALDB,  or the search for the 
 *     caldb file is successful. 
 *       
 *  2. dsGENERICERR :  this can happen only if myFile is CALDB.
 *     This error means that 
 *      1. if calProduct is not gaplookup, set myFile to NONE in main().
 *      2. if calProduct is gaplookup, 
 *            if 'flg' is INIT_NOT_OK, 
 *               set myFile to NONE in main().
 *            else 
 *               continue digging for the old degap file in main().
 *         
 * 
 * NOTEs:  
 *    1. calling sequence in main() should be : 
 *          getHdr -> init_caldb_var -> find_caldb_file ;
 *
 *    2. If input 'myFile' is not CALDB
 *          No change to myFile ;
 *          returns dsNOERR ;
 *          [ main() will use myFile directly. ]
 *
 *       else
 *          If All caldb function-calls are successfully, 
 *             Update myFile to the proper caldb file that's found ;
 *             return dsNOERR ;
 *             [ main() will use myFile directly. ]
 *          else 
 *             No change to myFile ; 
 *             return dsGENERICERR ; 
 *             [ main() will need to set myFile accordingly. ]
 *       endif
 *
 *    3. close 'myCaldb' in main(), not here !
 *--------------------------------------------------------------------*/
dsErrCode find_caldb_file(
                  char  *myFile,     /*U: see NOTEs #2 above. */
                  char  *myProduct,  /*I:*/
           HRC_CALDB4_P   hcp         /*U:*/
                    )
{
   int    nFile ;

  /*---------   input file is not CALDB   ----------*/
   if ( ds_strncmp_cis( myFile, "CALDB", 5) != 0 )   
      return dsNOERR ;

  /*----------  calInit was called before and it failed  --------*/
   if ( hcp->flg  == INIT_NOT_OK )
      return dsGENERICERR ;

   if ( (hcp->mySearch = calSetProduct(hcp->myCaldb, myProduct ) ) == NULL ) 
      return dsGENERICERR ;

   int Tmp = ds_map_hdr_to_caldb ( myFile, hcp->hdr, hcp->mySearch ) ;
   if ( Tmp != 0 )
      return dsGENERICERR ;

   if  ( ( nFile = calSearch ( hcp->mySearch) ) == 0 )
       return dsGENERICERR ;

   char *tmpStr;
   tmpStr = calGetFile(hcp->mySearch, 0); 
   if (tmpStr == NULL )
      return dsGENERICERR ;

  /*-------- successfully find the caldb file -----*/
   strcpy( myFile, tmpStr );
   calFree( tmpStr );

   if ( hcp->debug >= 2 )
      fprintf( stdout, "\nCALDB4 resolved to '%s'.\n", myFile );

   return dsNOERR ;

} /* end : find_caldb_file() */

/*-------------------------------------------------
 * If all files are not caldb, simply set the flag 
 * to INIT_NO_NEED and return.  Main() will not go 
 * thru any new caldb4 interface.
 *
 * use  hdrlib to get few keywords and
 * call calInit to set myCaldb and flg.
 *
 * Initialize mySearch to NULL.
 * Set detnam to hrc-s if the value is hrc-si.
 *--------------------------------------------------*/
HRC_CALDB4_P  init_caldb_var (
                INPUT_PARMS_P_T  inp_p      /* I:  */
                            )
{
    HRC_CALDB4_P  hcp = NULL ;
    hcp = (HRC_CALDB4_P) calloc( 1, sizeof(HRC_CALDB4_T)) ; 

    /*----  no need to have new caldb4 interface ----- */
    if ((ds_strncmp_cis( inp_p->degap_file, "caldb", 5) != 0 ) &&
        (ds_strncmp_cis( inp_p->gain_file, "caldb", 5) != 0 )  && 
        (ds_strncmp_cis( inp_p->adc_file, "caldb", 5) != 0 )  && 
        (ds_strncmp_cis( inp_p->hypfile, "caldb", 5) != 0 )  && 
        (ds_strncmp_cis( inp_p->tapfile, "caldb", 5) != 0 )  && 
        (ds_strncmp_cis( inp_p->ampflatfile, "caldb", 5) != 0 )  && 
        (ds_strncmp_cis( inp_p->ampsatfile, "caldb", 5) != 0 )  && 
        (ds_strncmp_cis( inp_p->ampsfcorfile, "caldb", 5) != 0 ) )
    {
       hcp->flg      = INIT_NOT_NEED ; 
       return hcp ;
    }

    hcp->myCaldb  = NULL ;
    hcp->mySearch = NULL ;
    hcp->flg      = INIT_NOT_OK ;
    hcp->debug    = inp_p->debug ;

    hcp->hdr = inp_p->caldb4_hdr;/*if use_obs=1, hdr from obs.par, else from 1st-evtin*/

    hdrGetKeyValue_c( hcp->hdr, "TELESCOP", &(hcp->telescop)) ;
    hdrGetKeyValue_c( hcp->hdr, "INSTRUME", &(hcp->instrume)) ;
    hdrGetKeyValue_c( hcp->hdr, "DETNAM", &(hcp->tmp_detnam )) ;

    if ( ds_strcmp_cis( hcp->tmp_detnam, "hrc-si" ) == 0 )
       strcpy( hcp->tmp_detnam, "hrc-s" ) ;

    if ( (hcp->myCaldb = calInit( hcp->telescop, hcp->instrume)) != NULL)
       hcp->flg  = INIT_OK ;

    return hcp ; 

} /*end: init_caldb_var() */
/*--------------------------------------*/
