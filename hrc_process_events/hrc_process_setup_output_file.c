/*                                                                
**  Copyright (C) 1997-2010  Smithsonian Astrophysical Observatory 
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
 
* FILE NAME: hrc_process_setup_output_file.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: The routine hrc_process_setup_output_file is responsible 
  for opening the output event file. It takes in a pointer to a structure 
  of EVENT_SETUP_T which contains all of the relevant evtlib data needed 
  from the input file to produce a 'new copy' for the output file. The 
  output file housekeeping information is then written to a second 
  EVENT_SETUP_T passed into the function. The routine does not return a 
  value. If any errors are detected they are stored in an error list and 
  retrieved latter via the error lib by the calling routine.

  JCC(8/23/01) - bugfix to avoid duplicate HISTORYs. 
  JCC(9/27/02) - write more CALDB filenames to header :
    HYPFILE  TAPRING  FLATFILE  SATFILE  AMPSFFIL
  JCC(3/2003)-replace old hdrlib with new one ( new putHdr ) 
             ( no Eventkey/PrimKey in ds_fix_evt_hdr.c )
JCC(4/2003)-add WIDTHRES key for tapRing
   (3/2004)-add RAND_SIZE key for rand_pix_size
(2/2005)-check block# before writing primary header; remove dmKernelGetCreate
(7/2005)-add removePath; write caldb files without paths;
    -add stkExpand for ASOLFILE; 
-Write keys for CALDB files even though they're NONE or blank (ADCCORF,BPIXFILE..).
 When it is 'none', convert it to NONE to prevent problem from merging L2.
(8/2005)-fixed stkExpand for '/path/a,/path/b' format
10/2009-dph/fap new gain files affect PI (see 'Notes on outCol PI')
1/2010-add hpeSetRang_s and hpe_set_ranges.
*H***********************************************************************/
 
/* hrc_process_events.h  includes delib.h */
#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#define HRC_PROCESS_EVENTS_H
#endif

#ifndef DS_HRC_CONFIG_H
#include "ds_hrc_config.h"
#define DS_HRC_CONFIG_H
#endif 

#ifndef HIST_HDR_H
#include "histlib.h"
#define HIST_HDR_H
#endif
 
void hrc_process_setup_output_file(
     EVENT_SETUP_P_T evtin_p,       /* I   - input file inprocession      */ 
     EVENT_SETUP_P_T evtout_p,      /* I/O - output file inprocession     */ 
     INPUT_PARMS_P_T inp_p,         /* I/O - input parameters/data        */ 
     ALIGNMENT_REC_P_T aln_p,       /* O   - alignment default values     */
     char*** colname,               /* O - output evt file col names  */
     dsErrList*      hpe_err_p)     /* O   - error list pointer           */
{
   char** out_names;                /* output evt file col names          */
   char * tmp = NULL ;
   boolean alloc_failure = FALSE;   /* indicates if an allocation failed  */ 
   int    yy=0 ;

   if (ds_clobber(evtout_p->file, (dsErrBool) inp_p->clobber, hpe_err_p)
       == dsNOERR)
   { 
      short num_cols;
 
      if ((evtout_p->extension = ds_create_event_file(evtout_p->file)) != NULL)
      {
        if ( dmBlockGetNo(evtout_p->extension) > 1 )
        {
          evtout_p->dataset =  dmBlockGetDataset(evtout_p->extension);
    /* evtout_p->primary=dmBlockCreate(evtout_p->dataset,"PRIMARY",dmIMAGE);*/
         evtout_p->primary=dmDatasetMoveToBlock(evtout_p->dataset, 1);
        }
        else
          evtout_p->primary = NULL ;
          
#ifdef NO_EVENTDEF
         evtout_p->num_cols = 25;

         /* allocate memory */
         evtout_p->mapping = (short*) calloc(evtout_p->num_cols, sizeof(short));
         evtout_p->desc = (dmDescriptor**) calloc(evtout_p->num_cols,
            sizeof(dmDescriptor*));
         evtout_p->types = (dmDataType*) calloc(evtout_p->num_cols,
            sizeof(dmDataType));
   
         if ((out_names = (char**) calloc(evtout_p->num_cols,
              sizeof(char*))) != NULL)
         {
            short rr;
   
            for (rr = 0; rr < evtout_p->num_cols; rr++)
            {
               if ((out_names[rr] = (char*)calloc(DS_SZ_COLUMN, sizeof(char)))
                    == NULL)
               {
                  alloc_failure = TRUE;
               }
            }
         }
   
         if ((evtout_p->desc == NULL) ||
             (evtout_p->mapping == NULL) ||
             (evtout_p->types == NULL) ||
             (out_names == NULL) ||
              alloc_failure)
         {
            /* set error status - memory allocation failed */
            dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Custom,
               "ERROR: Memory allocation failed setting up output event file.");
         }
         else
         {
            strcpy(out_names[0], "TIME");
            strcpy(out_names[1], "CRSU");
            strcpy(out_names[2], "CRSV");
            strcpy(out_names[3], "AMP_SF");
            strcpy(out_names[4], "AV1");
            strcpy(out_names[5], "AV2");
            strcpy(out_names[6], "AV3");
            strcpy(out_names[7], "AU1");
            strcpy(out_names[8], "AU2");
            strcpy(out_names[9], "AU3");
            strcpy(out_names[10], "RAWX");
            strcpy(out_names[11], "RAWY");
            strcpy(out_names[12], "CHIPX");
            strcpy(out_names[13], "CHIPY");
            strcpy(out_names[14], "TDETX"); 
            strcpy(out_names[15], "TDETY"); 
            strcpy(out_names[16], "DETX"); 
            strcpy(out_names[17], "DETY"); 
            strcpy(out_names[18], "X"); 
            strcpy(out_names[19], "Y"); 
            strcpy(out_names[20], "PHA");
            strcpy(out_names[21], "PI"); 
            strcpy(out_names[22], "SUMAMPS"); 
            strcpy(out_names[23], "CHIP_ID"); 
            strcpy(out_names[24], "STATUS");
   
            evtout_p->types[0] = dmDOUBLE;
            evtout_p->types[1] = dmSHORT;
            evtout_p->types[2] = dmSHORT;
            evtout_p->types[3] = dmSHORT;
            evtout_p->types[4] = dmSHORT;
            evtout_p->types[5] = dmSHORT;
            evtout_p->types[6] = dmSHORT;
            evtout_p->types[7] = dmSHORT;
            evtout_p->types[8] = dmSHORT;
            evtout_p->types[9] = dmSHORT;
            evtout_p->types[10] = dmSHORT;
            evtout_p->types[11] = dmSHORT;
            evtout_p->types[12] = dmSHORT;
            evtout_p->types[13] = dmSHORT;
            evtout_p->types[14] = dmLONG;
            evtout_p->types[15] = dmLONG;
            evtout_p->types[16] = dmDOUBLE;
            evtout_p->types[17] = dmDOUBLE;
            evtout_p->types[18] = dmLONG;
            evtout_p->types[19] = dmLONG;
            evtout_p->types[20] = dmSHORT;

          /*11/2009-outcol PI: maxPI=255||1023; SHORT dtype*/
          /*11/2009-   if ( inp_p->gainflag==OLD_SI_GAIN)  */
            evtout_p->types[21] = dmSHORT; /*outCol PI for new+old gainmap*/

            evtout_p->types[22] = dmSHORT;
            evtout_p->types[23] = dmSHORT;
            evtout_p->types[24] = dmLONG;
#else
         if (!ds_parse_eventdef(evtout_p->eventdef, &out_names, 
             &evtout_p->types, &num_cols )) 
         {
            evtout_p->num_cols = num_cols; 
            evtout_p->dim = (long*) calloc(evtout_p->num_cols, sizeof(long)); 
            evtout_p->mapping = (short*) calloc(evtout_p->num_cols, 
               sizeof(short));
            evtout_p->desc = (dmDescriptor**) calloc(evtout_p->num_cols,
               sizeof(dmDescriptor*));
         } 
         
         if ((evtout_p->desc == NULL) || (evtout_p->mapping == NULL) ||
             (evtout_p->dim == NULL) || (evtout_p->types == NULL) ||
             (out_names == NULL) || alloc_failure)
         {
            /* set error status - memory allocation failed */
            dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Custom,
               "ERROR: Memory allocation failed setting up output event file.");
         }
         else
         {
#endif 
            *colname = out_names ;   /*1/2010*/ 
            Header_Type *obs_info_p ;
           /* --- (9/2002) copy a full hdr. --- */
            yy=ds_copy_full_header(evtin_p->extension,evtout_p->extension,"hrc_process_events",1);
            if (yy != 0)     /* failed */
            {
               dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                       "ERROR: Failed to copy the full header\n");
            }

            /* write HDU keywords */
            /* JCC(2/2003) - moved to after the call of putHdr  */
   
            /*-----------------------------------------------------*/
            /* write out obs.par info (ie. obsfile) to the header  */
            /*-----------------------------------------------------*/
            if (inp_p->use_obs && inp_p->obs_info_p )
            {
               obs_info_p = inp_p->obs_info_p;     /* event hdr from obsfile */
   
               /* (9/2002) - fix some keys in event_hdr before writing it */
               /* (2/2003) - remove prim_info_p for new hdrlib in ds_fix_evt_hdr */
               ds_fix_evt_hdr(evtin_p->extension, obs_info_p) ;

               /*-- no HISTORY in obspar, always get them from 'infile' (ie. evtin_p) --*/
               /*   get_hist_lines (evtin_p->extension, &obs_info_p->histv,          */
               /*                   &obs_info_p->histnum, inp_p->debug);             */

               /* JCC(8/23/01)- set histnum=0 to avoid duplicate HISTORYs */
               /* obs_info_p->histnum = 0 ;   (2/2003) no histnum in new hdrlib */

/* JCC(2/2003) - replace with new hdrlib 'putHdr' */
               /* set_tool("hrc_process_events", USERTOOL);   - for CREATOR */
               /* put_event_hdr(evtout_p->extension, obs_info_p); */
               /* put_prim_hdr(evtout_p->primary, prim_info_p); */
               putHdr(evtout_p->extension, hdrDM_FILE, obs_info_p, EVENT_STS, "hrc_process_events");
           if ( evtout_p->primary != NULL ) 
               putHdr(evtout_p->primary, hdrDM_FILE, obs_info_p, PRIMARY_STS, "hrc_process_events");
            }
      
            /* write HDU keywords */
            /* JCC(2/2003) - moved to after the call of putHdr  */
            dmKeyWrite_c(evtout_p->extension, "CONTENT",
                         "EVT1", NULL, NULL);
            dmKeyWrite_c(evtout_p->extension, "HDUCLASS",
                         "OGIP", NULL, NULL);
            dmKeyWrite_c(evtout_p->extension, "HDUCLAS1",
                         "EVENTS", NULL, NULL);
            dmKeyWrite_c(evtout_p->extension, "HDUCLAS2",
                         "ALL", NULL, NULL);

            /* ---- write ASOLFILE to header ---- */
            if (inp_p->ASOLFILE)       /* is not null */
            {
               if (ds_strcmp_cis(inp_p->ASOLFILE, "NONE")==0) 
                  string_to_uppercase(inp_p->ASOLFILE, DS_SZ_PATHNAME);

               dmKeyWrite_c(evtout_p->extension, "ASOLFILE", inp_p->ASOLFILE,
                  NULL, "aspect offset files");
            }
            else
            {
               dmKeyWrite_c(evtout_p->extension, "ASOLFILE", "NONE",
                  NULL, "aspect offset files");
            }

            /* ---- write ADCfile to header ---- */
            if (ds_strcmp_cis(inp_p->adc_file, "NONE")==0) 
               string_to_uppercase(inp_p->adc_file, DS_SZ_PATHNAME);

               removePath(inp_p->adc_file, &tmp );
               dmKeyWrite_c(evtout_p->extension, "ADCCORF", tmp,
                  NULL, "ADC correction file used");

            /* ----- write bad pixel file to header ----- */ 
            if (ds_strcmp_cis(inp_p->badpixfile, "NONE")==0) 
               string_to_uppercase(inp_p->badpixfile, DS_SZ_PATHNAME);

               removePath(inp_p->badpixfile, &tmp) ;
               dmKeyWrite_c(evtout_p->extension, "BPIXFILE", tmp, 
                            NULL, "Bad pixel file used");

            /* ----- write gain file to header ----- */ 
            if (ds_strcmp_cis(inp_p->gain_file, "NONE")==0)  
               string_to_uppercase(inp_p->gain_file, DS_SZ_PATHNAME);

               removePath(inp_p->gain_file, &tmp) ;
               dmKeyWrite_c(evtout_p->extension, "GAINCORF", tmp,
                            NULL, "gain correction file used");

            /* ----- write degap file to header ----- */ 
            if (ds_strcmp_cis(inp_p->degap_file, "NONE")==0)  
               string_to_uppercase(inp_p->degap_file, DS_SZ_PATHNAME);

            removePath(inp_p->degap_file, &tmp) ;
            dmKeyWrite_c(evtout_p->extension, "DEGAP", tmp, NULL, NULL);
            if (ds_strcmp_cis(inp_p->degap_file, "COEFF") == 0)
            {
               dmKeyWrite_d(evtout_p->extension, "DEGAP_U1", inp_p->cf[0][0], 
                    NULL, "u degap linear correction factor");
               dmKeyWrite_d(evtout_p->extension, "DEGAP_U2", inp_p->cf[0][1], 
                    NULL, "u degap quadratic correction factor");
               dmKeyWrite_d(evtout_p->extension, "DEGAP_V1", inp_p->cf[1][0], 
                    NULL, "v degap linear correction factor");
               dmKeyWrite_d(evtout_p->extension, "DEGAP_V2", inp_p->cf[1][1], 
                    NULL, "v degap quadratic correction factor");
            } 

            /* ----- write hyperbolic test file to header ----- */
            if (ds_strcmp_cis(inp_p->hypfile, "NONE")==0)  
               string_to_uppercase(inp_p->hypfile, DS_SZ_PATHNAME);

               removePath(inp_p->hypfile, &tmp) ;
               dmKeyWrite_c(evtout_p->extension, "HYPFILE", tmp,
                            NULL, "hyperbolic test file used");

            /* ----- write tapring file to header ----- */
            if (ds_strcmp_cis(inp_p->tapfile, "NONE")==0)  
               string_to_uppercase(inp_p->tapfile, DS_SZ_PATHNAME);

               removePath(inp_p->tapfile, &tmp) ;
               dmKeyWrite_c(evtout_p->extension, "TAPRING", tmp,
                            NULL, "tapring coeff file used");

            /* ----- write flatness test file to header ----- */
            if (ds_strcmp_cis(inp_p->ampflatfile, "NONE")==0)  
               string_to_uppercase(inp_p->ampflatfile, DS_SZ_PATHNAME);

               removePath(inp_p->ampflatfile, &tmp) ;
               dmKeyWrite_c(evtout_p->extension, "FLATFILE", tmp,
                            NULL, "flatness test file used");
             
            /* ----- write saturation test file to header ----- */
            if (ds_strcmp_cis(inp_p->ampsatfile, "NONE")==0)  
               string_to_uppercase(inp_p->ampsatfile, DS_SZ_PATHNAME);

               removePath(inp_p->ampsatfile, &tmp) ;
               dmKeyWrite_c(evtout_p->extension, "SATFILE", tmp,
                            NULL, "saturation test file used");
             
            /* ----- write amp_sf correction file to header ----- */
            if (ds_strcmp_cis(inp_p->ampsfcorfile, "NONE")==0)  
               string_to_uppercase(inp_p->ampsfcorfile, DS_SZ_PATHNAME);

               removePath(inp_p->ampsfcorfile, &tmp) ;
               dmKeyWrite_c(evtout_p->extension, "AMPSFFIL", tmp,
                            NULL, "AMF_SF Correction file used");
             
            /* --- write WIDTHRES for tapRing if it's found in obs or evt ---*/
            if (inp_p->obs_widthres != NEG_9999 )
               dmKeyWrite_s(evtout_p->extension,"WIDTHRES",inp_p->obs_widthres,NULL,"width threshold");
            else if (inp_p->evt_widthres != NEG_9999 )
               dmKeyWrite_s(evtout_p->extension,"WIDTHRES",inp_p->evt_widthres,NULL,"width threshold");

            /* --- write RAND_SKY with the value of rand_pix_size ---*/
            dmKeyWrite_d(evtout_p->extension, "RAND_SKY", inp_p->randpixsize,
                    NULL, "pixel randomization");

      
            /* set up mappings */
            if (parse_hrc_evt_columns(out_names, evtout_p->num_cols, 
                evtout_p->mapping))  
            { 
               short bad_cols = evtout_p->num_cols;
   
               while (bad_cols--)
               {
                  if (evtout_p->mapping[bad_cols] == HDET_UNKNOWN_FIELD)
                  {
                     /* unable to map all requested output columns */
                     dsErrAdd(hpe_err_p, dsHPEBADOUTCOLERR, Individual, Generic,
                        out_names[bad_cols]);
                  } 
               } 
            }
            else
            { 
               /* JCC(2/2003) - remove EventKey from hrc_setup_columns for new hdrlib */
               hrc_setup_columns(evtout_p, inp_p, out_names, hpe_err_p); 
   
               /* write history comments */
	       put_param_hist_info(evtout_p->extension, "hrc_process_events", 
				   PFFile, inp_p->debug);

	       ds_write_pixhist_in_dm( evtout_p->extension );
               /*put_history(evtout_p->extension, "hrc_process_events",
                           evtin_p->file, evtout_p->file,
                           inp_p->debug);*/
   
            } /* end if parse_evt_column */ 
   
            /* free(out_names);  1/2010*/

         } /* end : evtout->desc/mapping/dim/types/out_names/alloc_failure */
      }  /* end : ds_create_event_file != NULL  */
      else
      {
         /* unable to create output event file */
         dsErrAdd(hpe_err_p, dsCREATEFILEERR, Individual, Generic, 
                  evtout_p->file);
      }
   } /* end: ds_clobber */
   /* printf("colname[1]=%s,\n", (*colname)[1]) ; */
} /* END :  hrc_process_setup_output_file */

/*---------------------------------------------------
 (7/2005) - remove pathnames from input string
 --------------------------------------------------- */
void removePath( 
    char*  inStr,     /* I */
    char** outStr     /* O */
    )
{
   char *ptr = NULL  ;

   ptr = strrchr( inStr, '/' ); 

   (*outStr) = (char*)calloc(DS_SZ_PATHNAME, sizeof(char)) ;  /* 4096 */

   if ( ptr != NULL )      /* find paths */
      strcpy( *outStr,  ptr+1 )  ;   
   else 
      strcpy( *outStr,  inStr )  ;   

   strtok( (*outStr), "[") ;

} /* END : strip off path */

/*------------------------------------------------------
 (7/2005) Expand stack to its content without any paths; 
          If inStr is NULL, return a NULL outStr.
 ------------------------------------------------------ */
void stkExpand( char *  inStr, /*I: stack||single_file||NULL*/
        char ** outStr /*O: string with filenames on stack*/
            )
{
      char    * atFile = NULL ;
      char    ppp[DS_SZ_PATHNAME] ;
      int     kkk=0 ;
      Stack   myStk=NULL ;

      if ( inStr == NULL )
      {
         *outStr = NULL ;
         return ;
      }

      (*outStr) = (char*)calloc(DS_SZ_PATHNAME, sizeof(char)) ;  /*4096*/

      memset( ppp, 0, DS_SZ_PATHNAME) ;

      myStk = stk_build( inStr );

      while( myStk && (( atFile = stk_read_next( myStk ) ) != NULL ))
      {
          char * tmp =NULL ;

          /* remove path from atFile */
          removePath ( atFile, &tmp ) ;

          strncat( ppp, tmp, strlen(tmp) );
          strncat( ppp, ",", 1 );
          free( atFile );
          atFile  = NULL;

      }
      stk_close( myStk );

      /* remove the last substring ", " */
      if ( (kkk=strlen(ppp)) > 1 )
          strncpy( (*outStr), ppp, kkk-1 ) ;  

} /* END: stkExpand */

/*1/2010 - initial version */
void hpe_set_ranges( EVENT_SETUP_P_T evtout_p,  /* I */
                     INPUT_PARMS_P_T inp_p,     /* I */
                              char** colname    /* I */
                   )
{
   short rr ;
   for (rr = 0; rr < evtout_p->num_cols; rr++)
   {
      switch (evtout_p->mapping[rr])
      {
         case HDET_PI:            /* 5/2010 */
            if ( inp_p->gainflag == OLD_SI_GAIN )
            {
               hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr], 
                              colname[rr], DS_HRC_PI_TLMIN, 
                              HDET_MAX_PI_OLD);  /*0<=short:pi<=255*/
            }
            else
            {
               hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr], 
                              colname[rr], DS_HRC_PI_TLMIN, 
                              HDET_MAX_PI_NEW);  /*0<=short:pi<=1023*/
            }
         break ;
         case HDET_PHA:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                   colname[rr], DS_HRC_PHA_TLMIN, DS_HRC_PHA_TLMAX);
         break ;
         case HDET_CP_X:
            switch (inp_p->hrc_system)
            {
               case HRC_IMG_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                       colname[rr], DS_HRCI_CRSU_TLMIN, DS_HRCI_CRSU_TLMAX);
               break;
               case HRC_SPC_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                        colname[rr], DS_HRCS_CRSU_TLMIN, DS_HRCS_CRSU_TLMAX);
               break;
               case HRC_SPC_IMG_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                         colname[rr], DS_HRCS_CRSU_TLMIN, DS_HRCS_CRSU_TLMAX);
               break;
               case HSI_IMG_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                          colname[rr], DS_HRCI_CRSU_TLMIN, DS_HRCI_CRSU_TLMAX);
               break;
               default:
               break;
            }
         break ; /* end: case HDET_CP_X */
         case HDET_CP_Y:
            switch (inp_p->hrc_system)
            {
               case HRC_IMG_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                          colname[rr], DS_HRCI_CRSV_TLMIN, DS_HRCI_CRSV_TLMAX);
               break;
               case HRC_SPC_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                          colname[rr], DS_HRCS_CRSV_TLMIN, DS_HRCS_CRSV_TLMAX);
               break;
               case HRC_SPC_IMG_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                          colname[rr], DS_HRCS_CRSV_TLMIN, DS_HRCS_CRSV_TLMAX);
               break;
               case HSI_IMG_SYS:
                  hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                          colname[rr], DS_HRCI_CRSV_TLMIN, DS_HRCI_CRSV_TLMAX);
               break;
               default:
               break;
            }
         break; /* end: case HDET_CP_Y */
         case HDET_AX_1:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                      colname[rr], DS_HRC_AU1_TLMIN, DS_HRC_AU1_TLMAX);
         break;
         case HDET_AX_2:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                      colname[rr], DS_HRC_AU2_TLMIN, DS_HRC_AU2_TLMAX);
         break;
         case HDET_AX_3:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                      colname[rr], DS_HRC_AU3_TLMIN, DS_HRC_AU3_TLMAX);
         break;
         case HDET_AY_1:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                      colname[rr], DS_HRC_AV1_TLMIN, DS_HRC_AV1_TLMAX);
         break;
         case HDET_AY_2:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                      colname[rr], DS_HRC_AV2_TLMIN, DS_HRC_AV2_TLMAX);
         break;
         case HDET_AY_3:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                      colname[rr], DS_HRC_AV3_TLMIN, DS_HRC_AV3_TLMAX);
         break;
         case HDET_AMP_SF:
            hpeSetRange_s( evtout_p->extension, evtout_p->desc[rr],
                      colname[rr], DS_HRC_AMP_SF_TLMIN, DS_HRC_AMP_SF_TLMAX);
         break;
      } /* end: switch (mapping) */
   }  /* end: for rr ..*/
} /*end: hpe_set_ranges */

/* 1/2010 - initial version */
void hpeSetRange_s(dmBlock* bb, dmDescriptor* dd, char *n1, short v1, short v2)
{
   short min = v1 ;
   short max = v2 ;
   dmDescriptorSetRange_s( dd, min, max);

   dmDescriptor* ss = dmSubspaceColOpen( bb, n1) ;
   if (ss)
   {
      /* dmSubspaceColUpdate_s( ss, &min, &max, 1 ); */
      dmSubspaceColSet_s( ss, &min, &max, 1 );
   }
} /* end: hpeSetRang_s */


