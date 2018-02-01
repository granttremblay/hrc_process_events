/*                                                                
**  Copyright (C) 1996-2010  Smithsonian Astrophysical Observatory 
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

* FILE NAME: hrc_process_events.c

* DEVELOPEMENT: tools

* DESCRIPTION:

  hrc_process_events program.  This tool is an HRC Level 1 processing tool.
  Its purpose is to read an event list, and for each event, determine
  its (tiled detector, detector, and/or sky) coordinates.

  The event list may be either a file or a 'stack'. The program outputs a
  qpoe file with columns defined in the hrc_process_events parameter file.
  A debug log may be selected as an option. If so, the user may also  
  specify whether to write the log to an output file or direct it to 
  stdout.

  hrc_process_events calls the following subroutines- sum_phas_hrc,
  calculate_coords_hrc, calculate_pi_hrc, write_hrc_events, 
  load_event_data, parse_hrc_evt_columns, dependency_check_init, 
  output_coord_validate, allocate_degap_table, degap_table_load, 
  deallocate_degap_table, load_bad_pixel_files, check_for_bad_pixels, 
  ratio_checks_hrc, allocate_adc_table, adc_table_load,
  apply_adc_correction, deallocate_adc_table, dependency_check_hrc,
  process_warnings, and cleanup_bad_pixel_data. 

* NOTES:

  wmclaugh@cfa  March 25, 1996  First Version.

* REVISION HISTORY:
  JCC(5/1/00) : updated for tap corrections 
     - add new functions : open_tap_ring_file() and check_tap_ring() 
     - add tap_ring_defs.h
     - update hrc_process_events.h to add tapfile, amps_3RD_raw, 
       and amps_tap_flag for EVENT_REC_P_T. 
     - update load_event_data() to initialize amps_3RD_raw as A3
     - add 'tapfile' to load_input_parameters() and hrc_process_events.par 
     - update write_hrc_events() to write amps_3RD_raw[u:vaxes] instead 
       of amps[][3]
     - add notes to calc_coarse_coords() and sum_phas_hrc().
     - move the following 3 calls outside of the loop: 
       open_amp_saturation_file; open_evt_flatness_file; open_hyperbolic_file;
     - 'check_tap_ring()' sets amp3[][3RD] between 4095 and 0. 
  JCC(5/11/00)
     - display the badpixfile warning before processing events.  
     - Keep the amplitudes as 'double' for all corrections : 
         add "double amps_dd" and "short amps_sh";  
         amps_dd is for computation and amps_sh is for outfile.
           adc_corr_routines.c  - use amp_dd for computation
           adc_filter_routines.c  - use amp_dd for computation
           coordinate_transforms.c  - use amp_dd for computation
           load_event_data.c - initialize amps_dd as amps_sh
           hrc_process_events.c - display amps_sh and amps_dd
           sum_phas_hrc.c - get the total amplitude from amps_dd
           tap_ring_functions.c - use amp_dd for computation
           write_hrc_events.c - write the raw amplitudes (amps_sh)
  JCC(5/12/00)-Saturation coeffs are a function of amp_sf;:
               rewrite the function "open_amp_saturation_file"  ;
               update the function "check_amp_saturation" ;
               change sat_test_coeffs_p to a structure.
  JCC(7/18/00)-set pi=pha if no gain map is used (do_pi).
              -add 1 to gain_index. 
  JCC(6/25/01) - if (PI > 255), set to 255.  (evt_p->pi)
               - made changes to calculate_pi_hrc as well.
  JCC(9/2002)- add amp_sf corrections; 
   -add new parameters to inp_p  :
      inp_p->do_amp_sf_cor ( from hpe.par )
      inp_p->range_switch_level(from obs.par or evt1.fits; write to outfile)
      inp_p->get_range_switch_level ( TRUE/FALSE ; intermediate param ; )
      inp_p->AMPSFCOR ( output TRUE/FALSE )
   - range_switch_level is also defined in the struct 'ampsfcor_coeff_str'
     which represents the column of range_switch_level in amp_sf_cor caldb.
   -add 4 new functions to amp_sf_cor_functions.c :
      open_amp_sf_cor_file() ;
      sum_raw_amps() ;
      apply_amp_sf_cor() ;
      write_amp_sf_corr() ;
   -regardless do_amp_sf_cor, always write out the key rangelev that comes 
    from obsfile or evt infile; obsfile supersedes evt .
 JCC(2/2003)- replace old hdrlib with a new one ( remove inp_p->primhdr_p) 
JCC(4/2003)-new tapRing spec ; 
JCC(7/2003)-reset status bit before reapplying any correction including
            tring,hyperbolic,amp_saturation,evt_flat (bit 0-5,30,31)
     -no change to the keyword ampsfcor for amp_sf corr.
     -add a new function 'initial_status' to reset other status bits;
 JCC(8/2003)- for an empty file, dmTableGetRowNo returns 1 instead of 
    dmBADROW in new dmlib; call dmTableSetRow to fix the problem.
 JCC(2/3/2004)-replace dmBlockIntersectSubspace with dmBlockCopy
 (6/2004)-condition check on pixlib 
 (4/2005)-freeHdr(obs_info_p);
 (2/2005,6/2005)-change printf format
(7/2005)-new stkExpand; expand ASOL stk without paths (ASOLFILE).
(8/2005)-add eFile to display #3 degap's warning with evt filename.
(5/2007)-fix segv with lev2 infile.
(5/2007)-new caldb4 ( calClose ) . 
(1/2009)- check infile-stack & infile-access at the beginning of hpe.c;
        - add hpePrintErr;
        - add peter hrcS 3dim gain image (obsolete 10/2009)
10/2009- replace peter hrcS 3dim gain image  w/ dph new hrcS gain table.
10/2009- add fap new hrcI gain image
11/2009- change outCol PI dtype from long to short for new gain files.
Notes on outCol PI (use OLD_SI_GAIN to set different output range ):
    for old 2dim gain image :    maxPI=255    dtype=short
    for new hrcS gain table :    maxPI=1023   dtype=short !
    for new hrcI gain image :    maxPI=1023   dtype=short !

5/2010 - write out ASPTYPE when infile has 0 row.
*H***********************************************************************/

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

#ifndef L1_ASPECT_DEFS_H
#include "l1_aspect_defs.h"
#define L1_ASPECT_DEFS_H
#endif

#ifndef ADC_COOR_DEFS_H
#include "adc_corr_defs.h"
#define ADC_COOR_DEFS_H
#endif 

#ifndef DS_HRC_CONFIG_H
#include "ds_hrc_config.h"
#define DS_HRC_CONFIG_H
#endif

dsErrCode hrc_process_events(void)
{
    /* CL INPUT PARAMETERS */ 
    INPUT_PARMS_T    inpars;    /* I - hold parameters used for processing */
    INPUT_PARMS_P_T  inp_p = &inpars; /* I - pointer to input parameters   */

    /* ALIGNMENT/ASPECT FILE VARIABLES */ 
    char  *align_file_p;        /* pointer to alignment file/stack name    */
    ALIGNMENT_INFRA_T  aln_infra; /* structure holding alignment file info */
    ALIGNMENT_INFRA_P_T aln_hk_p = &aln_infra; /* pntr to alignment HK info*/
    ALIGNMENT_REC_T  aln_rec;   /* data struct for alignment file record   */
    ALIGNMENT_REC_P_T aln_p = &aln_rec; /* ptr to an alignment file record */ 

    char            *asp_file_p; /* pointer to aspect file/stack name      */
    char *ASOLFILE=NULL ;  /* output the asol filenames without paths */

    ASPECT_INFRA_T   a_infra;   /* structure holding aspect file info      */
    ASPECT_INFRA_P_T asp_hk_p = &a_infra; /* pointer to aspect HK info     */
    ASPECT_ENTRY_T   asp_rec;   /* data struct to hold aspect file record  */
    ASPECT_ENTRY_P_T asp_p = &asp_rec; /* pointer to a aspect file record  */

    /* INPUT "STACK" VARIABLES */ 
    Stack evtstack = NULL;      /* input event file(s) 'stack'             */
    int    num_evtfile = 0;     /* number of event files in 'stack'        */ 

    /* OUTPUT EVENT FILE VARIABLES */ 
    EVENT_SETUP_T evt_out;      /* output event file information           */
    EVENT_SETUP_P_T evtout_p = &evt_out; /* pointer to output event file   */
    boolean setup_outfile = TRUE; /* fatal error has occurred if FALSE     */

    /* BAD EVENT FILE VARIABLES */ 
    EVENT_SETUP_T evt_bout;     /* output bad event file information       */
    EVENT_SETUP_P_T evtbout_p = &evt_bout; /* pointer to bad event file    */
    boolean setup_badfile = TRUE; /* T= bad event file has not been created*/

    /* STATISTICS/LOGFILE VARIABLES */
    STATISTICS_T stat;          /* struct to keep track of event statistics*/
    STATISTICS_P_T stat_p = &stat;/* pointer to statistics structure       */
    long      bad_interval = 0; /* keep track of # consecutive bad times   */ 
    register int     debug;     /* register reference to debug param value */ 
    FILE            *log_ptr;   /* L - pointer to logfile                  */ 

    /* DEGAP/ADC CORRECTION VARIABLES */ 
    DEGAP_CONFIG_P_T dgp_p = NULL;    /* degap table structure             */
    ADC_CORR_P_T adc_x = NULL;  /* pointer to x axis adc correction table  */
    ADC_CORR_P_T adc_y = NULL;  /* pointer to y axis adc correction table  */

    /* GAIN CORRECTION VARIABLES */
    float*      gain_p = NULL;  /* for old 2dim gain map image */
    EVENT_REC_T evt;            /* struct containing data for an event     */
    EVENT_REC_P_T evt_p = &evt; /* pointer to the event structure          */
    INST_KEYWORDS_T inst;       /*  instrument keyword data structure      */
    INST_KEYWORDS_P_T inst_p = &inst; /* instrument keyword structure pntr */
    double    last_time = -1.0; /* to keep track of out of sequence events */
    BAD_PIX_A_T  hotpix_p = {NULL, NULL, NULL}; /* hot spot (bad pixel) lists*/
    char*     evtfile = NULL;   /* pointer to store input file name        */
    dsErrCode erR = dsNOERR;    /* return error status for this routine    */
    dsErrList* hpe_err_p = NULL; /* pointer to list of errors/warnings     */

    char** e_names;   /*eventdef's col names */
    char** b_names;   /*bad event's col names */

    /* ADC FILTERING VARIABLES */
    double *flat_test_coeffs_p = NULL;
    /* short *sat_test_coeffs_p = NULL; */
    SAT_TEST_P_T       sat_test_coeffs_p = NULL ;
    HYP_TEST_P_T       hyp_test_coeffs_p = NULL;
    TRING_COEFFS_P_T   tring_coeffs_p = NULL;

    /* amp_sf_correction   */ 
    AMPSFCOR_COEFF_P_T ampsfcor_coeff = NULL ;

    /* initialize error list */ 
    dsErrCreateList(&hpe_err_p);

    /* initialize data structures */
    memset(evtout_p, 0, sizeof(EVENT_SETUP_T));
    memset(aln_p, 0, sizeof(ALIGNMENT_REC_T)); 
    memset(aln_hk_p, 0, sizeof(ALIGNMENT_INFRA_T)); 
    memset(asp_p, 0, sizeof(ASPECT_ENTRY_T)); 
    memset(asp_hk_p, 0, sizeof(ASPECT_INFRA_T)); 
    memset(stat_p, 0, sizeof(STATISTICS_T));

    /* load input parameters from 'hrc_process_events.par' */ 
    /* (4/2003)-intialized variables for inp_p */
    load_input_parameters(inp_p, hpe_err_p); 
    debug = inp_p->debug; 

    /* open output logfile or redirect to STDOUT */ 
    log_ptr = hrc_process_setup_logfile(inp_p, hpe_err_p); 

   /*------------------------------------------------------------------
    * (1/2009)-check infile-stack & infile-access here. If fatal, stop!
    *------------------------------------------------------------------*/
    evtstack = stk_build(inp_p->stack_in);
    short fatalE = 0 ;
    if ((num_evtfile = stk_count(evtstack)) < 1)
    {
       /* fatal error - empty stack! */
       dsErrAdd(hpe_err_p, dsSTKEMPTYERR, Individual, Custom,
              "ERROR: Input event stack is empty- expected a stack or file.");
       fatalE = 1 ;
    }
    else
    {
       while (((evtfile = stk_read_next(evtstack)) != NULL) && ( fatalE == 0))
       {
          if (dmDatasetAccess(evtfile,"R") != dmTRUE)
          {
             /* fatal error - can't access infile! */
             dsErrAdd(hpe_err_p, dsOPENFILEFERR, Individual, Generic, evtfile);
             fatalE = 1 ;      /* break ; */
          }
       }
    }
    stk_close(evtstack);   /* close stk here! */
    if (fatalE == 1 ) 
    {
        erR = hpePrintErr( hpe_err_p, log_ptr, inp_p->debug);
        return (erR);
    }
   /*------------------------------------------------------------------*/

    /* read obs.par file */
    /* 8/2002 - add to get inp_p->range_switch_level from obsfile ;
     *          Key name is RANGE_SWITCH_LEVEL in obsfile ; RANGELEV in evt ;
     */
    hrc_process_read_obsfile(inp_p, hpe_err_p); 

    /* setup calibration files and reset ADC and PI checks*/
    /* 8/2002 - add to get inp_p->range_switch_level from 
     *       event1 file if (obsfile==none) or (obsfile!=none, but it
     *       doesn't contain the keyname 'range_switch_level') ;
     *       Key name is RANGE_SWITCH_LEVEL in obsfile and RANGELEV in evt ;
     */
    access_caldb( inp_p, hpe_err_p );

    inp_p->do_ADC = ((ds_strcmp_cis(inp_p->adc_file, "NONE") != 0) && 
		     (strcmp(inp_p->adc_file, "\0") != 0)); 

    inp_p->do_pi = ((ds_strcmp_cis(inp_p->gain_file, "NONE") != 0) && 
		    (strcmp(inp_p->gain_file, "\0") != 0)); 

    /* if file is 'NONE' need to pass a NULL pointer to stack library*/ 
    if ((ds_strcmp_cis(inp_p->align_file, "NONE") == 0) || 
        (strcmp(inp_p->align_file, "\0") == 0))
    {
       align_file_p = NULL;
    }
    else
    {
       align_file_p = inp_p->align_file; 
    } 
    aln_hk_p->stack = stk_build(align_file_p);
    aln_hk_p->stk_cnt = stk_count(aln_hk_p->stack); 
    aln_hk_p->file = NULL; 

    /* if file is 'NONE' need to pass a NULL pointer to stack library */
    if ((ds_strcmp_cis(inp_p->asp_file, "NONE") == 0) ||
        (strcmp(inp_p->asp_file, "\0") == 0))
    {
       asp_file_p = NULL;
    }
    else
    {
       asp_file_p = inp_p->asp_file;

       /* expand stk to filenames without paths */
       stkExpand( asp_file_p, &ASOLFILE ) ;
       inp_p->ASOLFILE = ASOLFILE ;
    }
    /* if (inp_p->ASOLFILE != NULL) printf("ASOLFILE=%s\n",inp_p->ASOLFILE ); */
    asp_hk_p->stack = stk_build(asp_file_p);
    asp_hk_p->stk_cnt = stk_count(asp_hk_p->stack);
    asp_hk_p->asp_file = NULL;
    asp_p->curr = 0; 
    asp_p->next = 1; 


    /* 5/2010 - save ASPTYPE key in the 1st asol file here;  */
    char  asp_type[DS_SZ_KEYWORD] ;   
    memset( asp_type, 0 , DS_SZ_KEYWORD);
    if (asp_hk_p->stk_cnt >= 1 )
    {
       char* asol_p = NULL;   
       if (( asol_p = stk_read_num(asp_hk_p->stack, 1)) !=NULL)  /*1st file*/
       {
          dmBlock *aspBlock = dmTableOpen( asol_p );
          if ( aspBlock != NULL  )
             /* get the header key ASPTYPE from the first asol1 file */
             dmKeyRead_c( aspBlock, ASP_TYPE_KEY, asp_type, DS_SZ_KEYWORD);
          dmTableClose( aspBlock );
       }
    }


    if (debug > DEBUG_LEVEL_0)
    {
       fprintf(log_ptr," stack count is %2d items ...\n", num_evtfile); 
    }

    if (hpe_err_p->contains_fatal == 0)
    {
       if ((inp_p->start = parse_coord_range(inp_p->start_coord)) == 
           HDET_INV_VAL)
       {
          dsErrAdd(hpe_err_p, dsHPEBADCOORDRNGERR, Individual, Generic);
       }
       else if ((inp_p->stop = parse_coord_range(inp_p->stop_coord)) == 
                HDET_INV_VAL)
       {
          dsErrAdd(hpe_err_p, dsHPEBADCOORDRNGERR, Individual, Custom,
             "ERROR: The coordinate transformation ending point specified in the .par file is invalid."); 
       }
       else if (map_start_column(inp_p->start, &inp_p->x_col, &inp_p->y_col))
       {
          dsErrAdd(hpe_err_p, dsHPEBADCOORDRNGERR, Individual, Custom,
             "ERROR: The coordinate transformation starting point must be either coarse, chip, or tdet.");  
       }
    }
    /********************************************************************
     * (8/2002) - perform amp_sf corrections
     ********************************************************************/
    if (inp_p->do_amp_sf_cor==TRUE)
    {
       if (inp_p->get_range_switch_level == TRUE)
       {
         /* depending on the value of inp_p->range_switch_level,
          * we'll find the proper 6 column data in ampsfcorfile,
          * and saved in ampsfcor_coeff.
          */
          open_amp_sf_cor_file ( inp_p , &ampsfcor_coeff, hpe_err_p ) ;
       }
       else
       {
         /* this should never happend here   */
          dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
         "WARNING:  range_switch_level is missing in obsfile or event infile.  skip amp_sf correction.\n");
       }
    }
    /********************************************************************
     * JCC(5/1/00) - open_tap_ring_file opens the tap ring coefficients 
     *   file and get all the coefficients.  If the file is not specified, 
     *   tring_coeffs_p will be returned as NULL. 
     ********************************************************************/
    open_tap_ring_file(inp_p, &tring_coeffs_p,  hpe_err_p);

    /*******************************************
     * open ADC filtering test files and get the coefficients.
     * JCC(5/3/00) - place these 3 calls below.
     * JCC(5/12/00)- change sat_test_coeffs_p to a structure 
     *******************************************/
    open_amp_saturation_file(inp_p, &sat_test_coeffs_p, hpe_err_p);
    open_evt_flatness_file(inp_p->ampflatfile, &flat_test_coeffs_p, hpe_err_p);
    open_hyperbolic_file(inp_p->hypfile, &hyp_test_coeffs_p, hpe_err_p);

    /********************************************************************
     * start going through stack of infile          
     ********************************************************************/
    evtstack = stk_build(inp_p->stack_in);
    while (((evtfile = stk_read_next(evtstack)) != NULL) && 
           (hpe_err_p->contains_fatal == 0))
    {
       EVENT_SETUP_T evt_in;   /* input event file information */
       EVENT_SETUP_P_T evtin_p = &evt_in;
       long        row_check = dmSUCCESS;
       char  *tmp=NULL ;

       evtin_p->file = evtfile; 

       /* open input file */
       hrc_process_setup_input_file(evtin_p, inp_p, stat_p, hpe_err_p);

       /* set up instrument specific info */
       hrc_process_set_instrume(inp_p, hpe_err_p);

       if (inp_p->processing == HRC_PROC_FLIGHT)
       {
          aln_hk_p->processing = ALN_FLIGHT_DATA;
       }
       else
       {
          aln_hk_p->processing = ALN_XRCF_DATA;
       }
 
       if ((setup_outfile) && (hpe_err_p->contains_fatal == 0)) 
       {
          evtout_p->file = inp_p->outfile; 
          evtout_p->eventdef = inp_p->outcols; 

          hrc_process_configure_pixlib(evtin_p, inp_p, inst_p, aln_p, 
                                       hpe_err_p);
          if (hpe_err_p->contains_fatal!=0) break;

          /*(10/2009)-set up gain file. return 3 values of gainflag */
          load_gain_image(inp_p->gain_file, inp_p, &gain_p, hpe_err_p); 

          /*(10/2009)outCol PI will depend on inp_p->gainflag */
          hrc_process_setup_output_file(evtin_p, evtout_p, inp_p,  
                                        aln_p, &e_names, hpe_err_p); 
 
          /* write RANGELEV and AMPSFCOR for amp_sf_correction */ 
          write_amp_sf_corr_key(evtout_p->extension, inp_p);

          /* write instrument parameters as output header keywords */
          write_instrume_params(evtout_p->extension, inst_p);

          /* setup degap tables */
          hpe_setup_degap_file(inp_p, &dgp_p, hpe_err_p);

          if (inp_p->do_ADC)
          {
             if (!allocate_adc_table(inp_p, &adc_x, &adc_y, hpe_err_p))
             {
                adc_table_load(inp_p, adc_x, adc_y, hpe_err_p);
             }
          }

          /* 1/2009:  load_gain_image(inp_p->gain_file,inp_p,&gain_p,hpe_err_p);*/

          if (hpe_err_p->contains_fatal == 0)
          {
             if (debug > DEBUG_LEVEL_1)
             {
                fprintf(log_ptr," number of output columns = %2d\n",
                        evtout_p->num_cols);
             } /* if (debug) */
 
             if (output_coord_validate(evtout_p->mapping, 
                 evtout_p->num_cols, inp_p))
             {
                dsErrAdd(hpe_err_p, dsHPEOUTCOLUMNERR,  Individual, Generic);
             }
          }

          /* set up hot pixel list- set do_raw flag if hot pixel list exists */
          if (load_bad_pixel_files(inp_p->badpixfile, &hotpix_p[0]) != 0)
          {
             inp_p->do_raw = TRUE; 
          } 
          else if (ds_strcmp_cis(inp_p->badpixfile,"NONE")!=0)
          {                                            /* badpix != NONE */  
             dsErrAdd(hpe_err_p, dsHBBEMPTYINFILEERR,  Individual, Custom,
         "WARNING: Unable to load the bad pixel file- continuing without it.");
             process_warnings(hpe_err_p, log_ptr, debug); /*JCC:add(5/11/00)*/
          } 

          if (debug > DEBUG_LEVEL_4)
          {
             BAD_PIX_P_T hotspot_p = hotpix_p[0];
           
             fprintf(log_ptr, " HOT SPOTS (BAD PIXELS)\n");  
             while (hotspot_p != NULL)
             {
                fprintf(log_ptr, 
                   "        (%4ld, %4ld) to (%4ld, %4ld)\n", 
                   hotspot_p->x[0], hotspot_p->y[0], 
                   hotspot_p->x[1], hotspot_p->y[1]); 
                hotspot_p = hotspot_p->next; 
             } 
             fprintf(log_ptr, "\n\n");  
          }

         /*---------------------------------------------------------
          * intersect subspace- keeps gti's if rerunning 
          *
          * (2/3/04)-use dmBlockCopy instead to create the 1st subspace
          *--------------------------------------------------------*/
         /*  dmBlockIntersectSubspace(evtin_p->extension,evtout_p->extension,evtout_p->extension); */
          dmBlockCopy(evtin_p->extension, evtout_p->extension, "SUBSPACE");
          setup_outfile = FALSE; 
       }
       else
       {
          dmBlockMergeSubspace(evtin_p->extension, evtout_p->extension,
                               evtout_p->extension);
       } /* end: if (setup_outfile..) */

       if ((debug > DEBUG_LEVEL_3) &&
           (hpe_err_p->contains_fatal == 0)) 
       {
          fprintf(log_ptr,"  open file : %s\n", evtin_p->file);
       }

       /* add eFile to display #3 degap's warning with evt filename */
    if (hpe_err_p->contains_fatal == 0)
    {
       removePath(evtfile, &tmp );  /* evtfile!=NULL; rm 'path+filter';*/
       if (dgp_p!=NULL)
          strcpy(dgp_p->eFile, tmp);
    }

       if (hpe_err_p->contains_fatal == 0)
       {
          boolean in_time_exists = FALSE;
          boolean out_time_exists = FALSE;
          short rr = 0;

          if (debug > DEBUG_LEVEL_4)
          {
             fprintf(log_ptr, "   TIME       CRU CRV  AU1  AU2  AU3  AV1");
             fprintf(log_ptr, "  AV2  AV3   SUM  ");
             fprintf(log_ptr, "\n"); 
          }

          /*******************************************************/
          /* only loop if 1 or more rows in the input event file */
          /*******************************************************/
          if (dmTableGetRowNo(evtin_p->extension) != dmBADROW)
{
          row_check = dmTableSetRow(evtin_p->extension, 1); /*(8/2003)*/
          while ((row_check != dmNOMOREROWS) &&      /* while(evt_next_row)*/
                 (hpe_err_p->contains_fatal == 0)) 
          {
             /* initialize event record structure */
             memset(evt_p, 0, sizeof(EVENT_REC_T));

             evt_p->time = inp_p->default_time; 

             /*************************************
              * load data into event record 
              *************************************/
             load_event_data(evtin_p, evt_p);
             initial_status(inp_p, evt_p);

             /*************************************
              * (8/2002) - amp_sf corrections
              *************************************/
             if (( inp_p->do_amp_sf_cor == TRUE) &&
                 ( inp_p->get_range_switch_level == TRUE) &&
                 ( inp_p->match_range_switch_level == TRUE)  )
             {
                apply_amp_sf_cor(evt_p, ampsfcor_coeff) ;
             }

             /*********************************************************
              * JCC(5/1/00) -
              *   check_tap_ring  computes the corrected A3 and stored in 
              *   evt_p->amps_dd[u:v_axes][3]. The new A3 will be used to 
              *   compute the fine coordinates and for other tests.
              *   
              *   The fine coordinates will be computed in the existing 
              *   function calculate_coords_hrc. 
              *********************************************************/
             if (tring_coeffs_p != NULL)
                 check_tap_ring(inp_p, evt_p, tring_coeffs_p );

             /* keep track of earliest and latest event times */
             if (evt_p->time < inp_p->evt_tstart)
             {
                inp_p->evt_tstart = evt_p->time;
             }
             else if (evt_p->time > inp_p->evt_tstop)
             {
                inp_p->evt_tstop = evt_p->time;
             }

             /* if HRC-i flight data extra bit should be removed */
             if ((evt_p->cp[HDET_PLANE_Y] >= 64) &&
                 (inp_p->hrc_system == HRC_IMG_SYS))
             {
                evt_p->cp[HDET_PLANE_Y] -= 64;
             }

             if (inp_p->do_ADC)
             {
                apply_adc_correction(adc_x, adc_y, inp_p, evt_p); 
             }

             /********************************************************
	      * perform ADC filtering tests, if ARDs were provided.
              * hyperbolic test :
              ********************************************************/
	     if (hyp_test_coeffs_p != NULL)
	     {
	        check_hyperbolic(evt_p, hyp_test_coeffs_p);
	     }

             /*---------------------
	      * saturation tests
              *--------------------*/
	     if (sat_test_coeffs_p != NULL)
	     {
	        check_amp_saturation(evt_p, sat_test_coeffs_p);
	     }

	     if (flat_test_coeffs_p != NULL)
	     {
	        check_evt_flatness(evt_p, *flat_test_coeffs_p);
	     }

             if (evt_p->time < last_time)
             {
                dsErrAdd(hpe_err_p, dsHPEEVENTSEQERR, Accumulation, Generic,
                   evtin_p->file);
                evt_p->status |= HDET_SEQUENCE_STS; 
                stat_p->sequence_err++; 
                bad_interval++; 
             }
             else
             {
                if ((debug > DEBUG_LEVEL_3) && (bad_interval != 0))
                {
                   fprintf(log_ptr, 
                      "%3ld bad events occurred between time %8f and %8f\n",
                      bad_interval, last_time, evt_p->time); 
                   bad_interval = 0; 
                } 
                last_time = evt_p->time; 
             } 

             /* add time offset to event time for alignment sequence */
             evt_p->time += inp_p->time_offset;

             if (alignment_update(aln_p, aln_hk_p, evt_p->time))
             {
                dsErrAdd(hpe_err_p, dsHPEALIGNMENTERR, Individual, Generic); 
             }  

             /* remove time offset added to event time for alignment sequence*/
             evt_p->time -= inp_p->time_offset;

             if (aspect_update(asp_p, asp_hk_p, evt_p->time))
             {
                dsErrAdd(hpe_err_p, dsHPEASPECTERR, Individual, Generic); 
             }

             /* update statistical file counts */
             stat_p->total_events_in++;

             /* set next-in-line status bit if needed */
             if (inp_p->next_in_line)
             {
                evt_p->status |= HDET_NEXT_IN_LINE_STS; 
             }

             /*********************************************************
              *   compute the fine coordinates
              *********************************************************/
             if (calculate_coords_hrc(evt_p, inp_p, stat_p, 
                                      &asp_p->entry[asp_p->next],
                                      dgp_p, asp_hk_p->asp_file_type,
				      hpe_err_p))
             {
                /* reject event */ 
                if (setup_badfile)
                {
                   evtbout_p->file = inp_p->badfile; 
                   evtbout_p->eventdef = inp_p->badoutcols;
                   hrc_process_setup_output_file(evtin_p, evtbout_p, 
                                          inp_p, aln_p, &b_names, hpe_err_p);
                   setup_badfile = FALSE; 
                } 

                /* update- add current event to bad event file */ 
                write_hrc_events(evtbout_p, evt_p, hpe_err_p); 
                dsErrAdd(hpe_err_p, dsHPEBADEVTFILEERR, Accumulation, Generic,
                         evtbout_p->file); 
             }      
             else 
             {
                if (debug > DEBUG_LEVEL_4)
                {
                   fprintf(log_ptr, 
                      "%9.4f   %3d %3d %4d %4d %4d %4d %4d %4d %5d\n",
                      evt_p->time, evt_p->cp[HDET_PLANE_X], 
                      evt_p->cp[HDET_PLANE_Y],
                      evt_p->amps_sh[HDET_PLANE_X][HDET_1ST_AMP],
                      evt_p->amps_sh[HDET_PLANE_X][HDET_2ND_AMP],
                      evt_p->amps_sh[HDET_PLANE_X][HDET_3RD_AMP],
                      evt_p->amps_sh[HDET_PLANE_Y][HDET_1ST_AMP],
                      evt_p->amps_sh[HDET_PLANE_Y][HDET_2ND_AMP],
                      evt_p->amps_sh[HDET_PLANE_Y][HDET_3RD_AMP], 
                      evt_p->sum_amps);
                }

                if (debug > DEBUG_LEVEL_4 )
                {
                   printf(
                      "%9.4f %6ld %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
                      evt_p->time, row_check, 
                      evt_p->amps_dd[HDET_PLANE_X][HDET_1ST_AMP],
                      evt_p->amps_dd[HDET_PLANE_X][HDET_2ND_AMP],
                      evt_p->amps_dd[HDET_PLANE_X][HDET_3RD_AMP],
                      evt_p->amps_dd[HDET_PLANE_Y][HDET_1ST_AMP],
                      evt_p->amps_dd[HDET_PLANE_Y][HDET_2ND_AMP],
                      evt_p->amps_dd[HDET_PLANE_Y][HDET_3RD_AMP],
                      evt_p->amp_tot[HDET_PLANE_X]);
                }

                /* use the gain map to calculate pi , or set pi=pha*/
                if (inpars.do_pi)
                {
                   calculate_pi_hrc(gain_p, inp_p, evt_p);
                }  
                else                                /* 7/18/00 */
                {
                   evt_p->pi = evt_p->pha ;

                  /*10/2009- see 'Notes on outCol PI' */
                   long HDET_MAX_PI_VALUE =  HDET_MAX_PI_OLD ;   /* 255 */
                   if (inp_p->gainflag != OLD_SI_GAIN )
                      HDET_MAX_PI_VALUE = HDET_MAX_PI_NEW ;      /* 1023*/

                   /* (6/25/01) - add */
                   if (evt_p->pi > HDET_MAX_PI_VALUE)     /*255 or 1023*/
                   {
                      evt_p->status |= HDET_PI_VALUE_STS; /*flag if pi > max */
                      evt_p->pi = HDET_MAX_PI_VALUE ; 
                   }
                   /* end: (6/25/01) */
                }

                /* set status bits of hot spots (bad pixels) */
                check_for_bad_pixels(hotpix_p, evt_p); 

                /* write data to output event file */
                write_hrc_events(evtout_p, evt_p, hpe_err_p);

                /* update statistical file counts */
                stat_p->total_events_out++;
             } /* end:  if (calculate_coords == TRUE) */ 
             row_check = dmTableNextRow(evtin_p->extension);
          } /* end:  while (evt_next_row)  */ 
} /* end : if ( dmTableGetRowNo != dmBADROW ) */
          /*******************************************************
           * end: loop thru all records of one infile
           *******************************************************/

          while ((rr< evtin_p->num_cols) && (!in_time_exists))
          {
             in_time_exists = (evtin_p->mapping[rr++] == HDET_TIME);
          }
          rr = 0;
          while ((rr < evtout_p->num_cols) && (!out_time_exists))
          {
             out_time_exists = (evtout_p->mapping[rr++] == HDET_TIME);
          }
 
          if (inp_p->end_time == inp_p->default_time)
          {
             inp_p->end_time += TIME_INT_EPSILON;
          }
       }  /* end: if (hpe_err_p->contains_fatal == 0) */ 

       if (process_warnings(hpe_err_p, log_ptr, debug))
       {
          stat_p->num_bad_files++;    
       } 
  
       if ((evtin_p->file != NULL) && (debug > DEBUG_LEVEL_0))
       {
          fprintf(log_ptr,"  close file: %s\n", evtin_p->file);
       }
       hrc_process_evt_file_cleanup(evtin_p);

       if (evtfile != NULL)
       {
          free(evtfile); 
          evtfile = NULL;
       }
    } /* end while (evtin_p->file) */ 
    /*******************************************************/
    /*  end: loop of the stack of infiles                  */
    /*******************************************************/
                          
    /* error if there are no good files in stack */
    if ((stat_p->num_files_in == stat_p->num_bad_files) && 
        (stat_p->num_files_in > 0))
    {
       dsErrAdd(hpe_err_p, dsHPEALLBADFILESERR, Individual, Generic); 
    }

    /* clean up degap tables */
    l1h_deallocate_degap_table(&dgp_p, hpe_err_p);

    /* free up memory from adc correction tables */ 
    deallocate_adc_table(&adc_x, &adc_y);

    /* free up memory allocated for hot pixel list */
    cleanup_bad_pixel_data(hotpix_p); 

    /* free memory for alignment/aspect files */
    close_alignment_file(aln_hk_p); 
    close_aspect_file(asp_hk_p); 
    if (asp_hk_p->asp_type[0] != '\0')
    {
       dmKeyWrite_c(evtout_p->extension, ASP_TYPE_KEY,
                    asp_hk_p->asp_type, NULL, NULL);
    }
    else
    {
     /*5/2010 both asp_type & asp_hk_p->asp_type come from ASPTYPE in asol.*/
      if ( asp_type[0] != '\0')
         dmKeyWrite_c(evtout_p->extension, ASP_TYPE_KEY, asp_type, NULL, NULL);
    }

    /* free up memory for old gain map */
    if ((gain_p != NULL) && (inp_p->gainflag==OLD_SI_GAIN)) 
    {
      free(gain_p);
      gain_p = NULL; 
    }
    /* 10/2009 - free up memory for new hrcS gain table*/
    if ( inp_p->gainflag==NEW_S_GAIN ) 
    {
       free(inp_p->gainmapVal);
       free(inp_p->tgainVal);
       free(inp_p->rawxgridVal);
       free(inp_p->rawygridVal);
       free(inp_p->timegridVal);
       free(inp_p->obs_tgain);
       free(inp_p->G_2nd);
    }

    /* free memory for tap ring corrections */
    if (tring_coeffs_p != NULL)
       free(tring_coeffs_p) ;

    /* free memory for ADC filtering tests */
    if (hyp_test_coeffs_p != NULL)
    {
       free(hyp_test_coeffs_p);
    }

    if (flat_test_coeffs_p != NULL)
    {
       free(flat_test_coeffs_p);
    }

    if (sat_test_coeffs_p != NULL)
    {
       free(sat_test_coeffs_p);
    }

    /* free up memory allocated for old_system */
    if (inp_p->old_sys != NULL)
    {
       free(inp_p->old_sys);
    } 

    /* close output file */
    hpe_set_ranges( evtout_p, inp_p, e_names ); 
    hrc_process_evt_file_cleanup(evtout_p); 

    /* if a bad event file was generated- do housekeeping */
    if (!setup_badfile)
    {
       hrc_process_evt_file_cleanup(evtbout_p);
    } 
    
    /* free up memory allocated for stack */
    stk_close(evtstack);

    /* free up memory allocated for pixlib */ 
    if (inp_p->pix_init)
    {
       pix_close_pixlib();  
       inp_p->pix_init = FALSE; 
       if (inp_p->obs_info_p) 
       {
          /* free any memory used for events extension keywords */
          freeHdr(inp_p->obs_info_p);
          inp_p->obs_info_p = NULL; 
       } 
    } 

    /* verify time range */
    if (hpe_err_p->contains_fatal == 0)
    {
       hrc_process_time_check(inp_p, log_ptr, hpe_err_p);
    }

    if (debug > DEBUG_LEVEL_2)
    {
       fprintf(log_ptr,"\n ============ STATISTICS ===========\n"); 
       fprintf(log_ptr,"FILES   total in (%4ld) - bad (%2ld) = %3ld used\n", 
          stat_p->num_files_in, stat_p->num_bad_files,
          stat_p->num_files_in - stat_p->num_bad_files); 
       fprintf(log_ptr,"EVENTS  total in (%4ld)   total out = %3ld\n", 
          stat_p->total_events_in, stat_p->total_events_out);  
       fprintf(log_ptr,"SPILLOVERS  neg (%3ld)   pos (%3ld)    total = %4ld\n", 
          stat_p->fixed_mfinpos, stat_p->fixed_pfinpos, 
          (stat_p->fixed_mfinpos + stat_p->fixed_pfinpos) ); 
       fprintf(log_ptr, "OUT of SEQUENCE events = %3ld\n",
          stat_p->sequence_err); 
    }

    erR = hpePrintErr( hpe_err_p, log_ptr, inp_p->debug);   /* 1/2009 */

    /* close up log file */
    if ((debug > DEBUG_LEVEL_0) && (log_ptr != NULL))
    {
       fclose(log_ptr); 
    } 

   if (inp_p->hcp->flg == INIT_OK ) 
      calClose(inp_p->hcp->myCaldb);

    return (erR);

}
