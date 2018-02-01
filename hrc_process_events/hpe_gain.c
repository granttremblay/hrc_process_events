/*                                                                
**  Copyright (C) 2009,2012  Smithsonian Astrophysical Observatory 
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


/*----------------------------------------------------------
 * This file contains functions used for dph hrcS gain table.
 *
 * JCC(10/2009) - initial version 
 * JCC(8/2012) - make TIMEGRID_LEN, RAWX_LEN dynamic for hrcS t_gain_map.
 *     ( Note: the old 'fixed' values were TIMEGRID_LEN=18, RAWX_LEN=48 )
 *----------------------------------------------------------*/

#include "hrc_process_events.h"
#include "ds_hrc_config.h"
#include "dsnan.h"

/*----------------------------------------------------------
 * 1. Use evt rawx/rawy to get the index of gainmap column
 *    from dph new hrcS gain table. 
 * 2. Compute the output PIn value.
 *---------------------------------------------------------- */
void S_new_gain_index_pi(
            INPUT_PARMS_P_T inp_p,    /*I: */
                EVENT_REC_T *evt_p    /*U: */
                   )
{
   short Fnd ;

  /* dph spec eq(7);  XJ ; 
   * Similar to evt_p->gain_index[0] used by old 2dim gain img*/

   int RAWX_LEN  =  ( int ) inp_p->rawxgridSize ;   /* 8/2012 */

   long xj = grid_idx( RAWX_LEN, inp_p->rawxgridVal, 
                       evt_p->rawpos[HDET_PLANE_X], &Fnd );

  /* dph spec eq(6) ;  YI ; 
   * Similar to evt_p->gain_index[0] used by old 2dim gain img*/
   long yi = grid_idx( RAWY_LEN, inp_p->rawygridVal, 
                       evt_p->rawpos[HDET_PLANE_Y], &Fnd );

  /* dph spec eq(8) ;  compute PIn */

   /* 8/2012: pass RAWX_LEN */
   evt_p->pi_double = inp_p->gainmapVal[gmIdx(0,xj,yi,RAWX_LEN)] + 
                        inp_p->G_2nd[g2Idx(xj,yi,RAWX_LEN)]*evt_p->DDn ;

  /*------------------------------------------------------------------------
   *Note: due to do_pi and outCol statusBit setup, checking on the range of 
   *      pi_double won't be done here! Instead, do it in calculate_pi_hrc()
   *------------------------------------------------------------------------*/

} /* end: S_new_gain_index_pi() */ 

/*------------------------------------------------------------------------
 * Use evt_mjd_obs to find the matched timegrid from dph hrcS gain table:
 *   (timeGrid=timegridVal[ww]) <= evt_mjd_obs < timegridVal[ww+1]
 * Return  timeGrid ;
 *------------------------------------------------------------------------*/
void  find_timegrid_val (
        INPUT_PARMS_P_T  inp_p,     /*I*/
          long           *ww,       /*O*/
          double         *timeGrid, /*O*/
          short          *found     /*O*/
              )
{

int TIMEGRID_LEN = (int) inp_p->timegridSize ; 

   *ww = grid_idx(TIMEGRID_LEN, inp_p->timegridVal, inp_p->evt_mjd_obs, found);/*0-17*/

   timeGrid[0] = inp_p->timegridVal[ *ww] ;

   if ( *found == 1 )
      timeGrid[1] = inp_p->timegridVal[ *ww +1 ] ;
   else 
      timeGrid[1] = timeGrid[0] ;

   return ;
}

/* -----------------------------------------------------------------------------
 * Pass evt rawx,rawy,mjd_obs and find the associated array index for  
 * the columns of rawxgrid/rawygrid/timegrid in dph hrcS gain table.
 *
 * rawxgrid/rawygrid/timegrid in dph hrcS gain table must be in ascending order.
 *
 * if evt rawx/rawy/mjd_obs is less than the min of rawxgrid/rawygrid/timegrid, 
 *        return idx = first_index = 0  ( set  found=0 )
 *
 * if evt rawx/rawy/mjd_obs is greater than or equal to the max of 
 *  rawxgrid/rawygrid/timegrid, 
 *        return idx = last_index       ( set  found=0 )
 *        ( last_index of rawxgrid = 47 )
 *        ( last_index of rawygrid = 575 )
 *        ( last_index of timegrid = 17 )
 *
 * otherwise, set found=1  and return idx : 
 *          rawxgrid[ idx ] <= evt_rawx    < rawxgrid[ idx+1 ]
 *          rawygrid[ idx ] <= evt_rawy    < rawygrid[ idx+1 ]
 *          timegrid[ idx ] <= evt_mjd_obs < timegrid[ idx+1 ]
 *
 *
 * 1. dim_Grid of rawxgrid/rawygrid/timegrid are ( 48, 576, 18 )
 * 2. gridCol can be the columns of rawxgrid/rawygrid/timegrid 
 * 3. return idx [ranged from 0-17(rawxgrid); 0-575(rawxgrid); 0-17(timegrid);]
 * ---------------------------------------------------------------------------*/
long grid_idx (
     long   dim_Grid,  /*I: dim of gridCol; 48/576/18 */
     double *gridCol,  /*I: array of grid columns in new hrcS gain table */
     double evtVal,    /*I: each evt's rawx,rawy;  or evt's key mjd_obs; */
     short  *found     /*O: 0=not-found; 1=found; */
            )
{
   long idx=0 ;    /* return array-index */

   *found = 0 ;                               /* not found */
   if ( dim_Grid == 1 )
      return 0 ;                              /* idx = 0 */

   if ( hpe_gt( gridCol[0], evtVal ) == 1 )   /* gridCol[0] > evtVal*/
      return 0 ;                              /* idx = 0 */

   if (hpe_ge(evtVal,gridCol[dim_Grid-1])==1) /* evtVal>=gridCol[dim_Grid-1])*/
      return (dim_Grid-1) ;                   /* idx = dim_Grid-1 */

   long   beg, end, mid ;    /* 'begin,end,mid' index */
   beg = 0 ;                 /* beg_rawx=0;  beg_rawy=0;   beg_timegrid=0 */ 
   end = dim_Grid - 1 ;      /* end_rawx=47; end_rawy=575; end_timegrid=17 */

   mid = get_mid ( beg, end ) ;

  /*------------------------------------------------------------------
   * search the proper index (NN) among gridCol to meet the condition below:
   *
   *                gridCol[ NN ]  <= evtVal  < gridCol[ NN+1 ] 
   * hpeDPHgain:   RAWXGRID[ NN ] <= evt_RAWX  < RAWXGRID[ NN+1 ] 
   * hpeDPHgain:   RAWYGRID[ NN ] <= evt_RAWY  < RAWYGRID[ NN+1 ] 
   * hpeDPHgain:   TIMEGRID[ NN ] <= evt_MJD_OBS  < TIMEGRID[ NN+1 ] 
   *
   * gridCol array is in ascending order ;
   *------------------------------------------------------------------*/
   while ((mid >= beg ) && (mid <= end) && (*found == 0 ))
   {
      if ( evtVal >= gridCol[ mid ] )
         beg = mid ;    
      else
         end = mid ; 

      mid = get_mid ( beg, end ); 

      if ( ( mid == -1 ) || ( beg == end )) 
      {
         if ( evtVal >=  gridCol[ end ] )
         {
            *found = 1 ;
            mid = end ;
         }
         else 
         {
            if ( evtVal >=  gridCol[ beg ] )
            { 
               *found = 1 ;
               mid = beg ;
            }
            else
            {
               *found = 0 ; 
               mid = -1 ;
               break ;
            }
         }
      } /* end : if mid== -1  */
   } /* end :  while 'mid' */

   if ( *found==0 )
   {
      idx = 0 ;         /* should never happen */
   }
   else if (*found==1 )
   {
        idx = mid ;      /* gridCol[ mid] <= evt_raw < gridCol[ mid +1 ] */
   }
   return idx ;
} /* end : grid_idx() */

/*-----------------------------------------------------------------------------
* Compute DDn coefficient from sum_amps and amp_sf.  ( ie. eq.5 in dph spec )
*
* ushort evt_p->sum_amps = raw_au1+raw_au2+au3_tap_corr+raw_av1+raw_av2+av3_tap_corr
* short  evt_p->amp_sf   = 'amp_sf' after amp_sf_corr applied
* double evt_p->DDn = evt_p->sum_amps * ( 2**(evt_p->amp_sf - 8) ) ;
* ---------------------------------------------------------------------------*/
void calc_DDn (
         EVENT_REC_T     *evt_p      /*I:sum_amps & amp_sf; O: DDn; */
               )
{
  /* dph spec eq.5 */
  evt_p->DDn = ldexp ( evt_p->sum_amps, evt_p->amp_sf - 8) ;

} /* end: calc_DDn() */

/* --------------------------------------------------------------
 * This function is used when gain_file is not NONE or blank.
 *
 * If gain_file is a table type, set flag to NEW_S_GAIN 
 * Else, set flag to OLD_SI_GAIN for now. load_gain_image will
 * check on SAMPNORM and may change flag to NEW_I_GAIN.
 * -------------------------------------------------------------*/
void set_gainflag(
      INPUT_PARMS_P_T inp_p,       /* I/O - gain_file, gainflag */
      dsErrList*      err_p        /* U */
                 )
{
   dmBlock *myB = NULL; 
   if (NULL == (myB = dmBlockOpen(NULL, inp_p->gain_file))) 
   {
      dsErrAdd (err_p, dsOPENFILEFERR, Individual, Custom, "ERROR: can't open gain file, %s\n", inp_p->gain_file );
      return ;
   }

   dmBlockType  myT ;
   myT = dmBlockGetType ( myB ) ;

   if ( myT == dmIMAGE )    /*2 dim gain image*/
   {
      inp_p->gainflag = OLD_SI_GAIN; /*can't distinguish it w/ NEW_I_GAIN here*/
   }
   else
   {
      inp_p->gainflag = NEW_S_GAIN; /* new hrcS gain table*/

     /*--- error out if mjd_obs not found in evt file */
      if (inp_p->mjd_obs_warn == 1 )
      {
         dsErrAdd (err_p, dsFINDKEYWORDFERR, Individual, Custom, "ERROR : mjd_obs keyword is missing in infile and can't be computed neither.\n");

         dmBlockClose( myB ) ;
         return ;
      }
      else if (inp_p->mjd_obs_warn == 2)  /*compute it successfully*/
         err_msg("WARNING: mjd_obs keyword is missing in infile. value from computation.\n");
   }
   dmBlockClose( myB ) ;
   return ;
      
} /* end: set_gainflag  */

/*----------------------------------------
 * read a column from new hrcS gain table 
 *----------------------------------------*/
void read_gainTab_col (
          dmBlock*  srcBlock,     /* I: */
          char*     colName,      /* I: col name*/
          double**  colVal,       /* O: array value*/
          long*     colSize       /* O: array size*/
                 )
{
    dmDescriptor* nameDesc ;

    if (srcBlock)
    {
       nameDesc = dmTableOpenColumn (srcBlock, colName);
       *colSize = dmGetArraySize( nameDesc ) ;
    }
    *colVal =  (double *)calloc( *colSize, sizeof(double))  ;
    dmGetArray_d( nameDesc, *colVal,  *colSize );

    long ii;
    for (ii=0; ii< (*colSize); ii++)
    {
       if ds_dNAN(((*colVal)[ii]))   /* it's nan */
       {
          (*colVal)[ii] = 0.0 ;
       }
    }

    return ;
} /* read_gainTab_col() */


/* ---------------------------------------------------
 * Read columns from dph new hrcS gain table :
 *    GAINMAP[2,48,576] = Real8
 *    TGAIN[576,18]     = Real8
 *    RAWXGRID[48]      = Real8
 *    RAWYGRID[576]     = Real8
 *    TIMEGRID[18] MJD  = Real8
 * ---------------------------------------------------*/
void load_S_new_gain_table(
            INPUT_PARMS_P_T inp_p   /* U */
                    )
{
   if ( inp_p->gainflag != NEW_S_GAIN )   
      return ;

   dmBlock* srcBlock = NULL ;
   srcBlock = dmTableOpen( inp_p->gain_file ) ;

  /* -----------------------------------------------------
   *          BEGIN: read new hrcS gain table 
   * -----------------------------------------------------*/
  /* --- column GAINMAP[ 2, 48, 576 ] in random order --- */
   inp_p->gainmapSize = 0 ;
   inp_p->gainmapVal = NULL ;
   read_gainTab_col( srcBlock, "GAINMAP", &inp_p->gainmapVal, &inp_p->gainmapSize);

  /* --- column TGAIN[ yi=576, ww=18 ] in decending order for a fixed yi --- */
   inp_p->tgainSize = 0 ;   /* 1dim array size */
   inp_p->tgainVal = NULL ;  /* 1dim array elements */
   read_gainTab_col( srcBlock, "TGAIN", &inp_p->tgainVal, &inp_p->tgainSize );

  /* --- column RAWXGRID[ 48 ] in ascending order --- */
   inp_p->rawxgridSize = 0 ;   /* 1dim array size */
   inp_p->rawxgridVal = NULL ;  /* 1dim array elements */
   read_gainTab_col( srcBlock, "RAWXGRID", &inp_p->rawxgridVal, &inp_p->rawxgridSize );

  /* --- column RAWYGRID[ 576 ] in ascending order --- */
   inp_p->rawygridSize = 0 ;   /* 1dim array size */
   inp_p->rawygridVal = NULL ;  /* 1dim array elements */
   read_gainTab_col( srcBlock, "RAWYGRID", &inp_p->rawygridVal, &inp_p->rawygridSize );

  /* --- column TIMEGRID[ ww=18 ] in ascending order --- */
   inp_p->timegridSize = 0 ;   /* 1dim array size */
   inp_p->timegridVal = NULL ;  /* 1dim array elements */
   read_gainTab_col( srcBlock, "TIMEGRID", &inp_p->timegridVal, &inp_p->timegridSize );

  /* -----------------------------------------------------
   *          END: read new hrcS gain table 
   * -----------------------------------------------------*/

  /* -----------------------------------------------------------
   * dph spec step3:  interpolate tgain to observed date. eg.(3)
   *                  (ie. compute obs_tgain[yi=0->575] )
   * Expect the values of obs_tgain[0->575] not to be in order.
   * -----------------------------------------------------------*/
   short fnd ;
   long  ww;              /* 0-17; for timegrid column in new hrcS gain table*/
   double timeGrid[2] ;   /* TIMEGRID[ww],[ww+1] in dph spec eq 3.
                           * Function of evt_mjd_obs  */

  /* use evt_mjd_obs to find proper timegrid in new hrcS gain table */
   find_timegrid_val( inp_p, &ww, timeGrid, &fnd ) ;  

   double a1 = 0.0 ; 
   double a2 = 1.0 ;       /* initial to 1 */
   a1 = inp_p->evt_mjd_obs - timeGrid[0] ; 
   if (fnd==1)
      a2 = timeGrid[1] - timeGrid[0] ;   /*if fnd==1, expect a2>0 */

   inp_p->obs_tgain=(double *)calloc(RAWY_LEN,sizeof(double));/*Tobs in dph eq.3*/
   long yi ;
   for (yi=0; yi<RAWY_LEN; yi++)    /* 0->575 */
   {
      double t_yi_ww0;                      /* tgain[yi,ww0] in dph spec eq.3*/
      double t_yi_ww1;                      /* tgain[yi,ww1] in dph spec eq.3*/

      t_yi_ww0  = inp_p->tgainVal[tgIdx(yi,ww)];         /* tgain[yi,ww0] */
      inp_p->obs_tgain[yi] = t_yi_ww0 ;                     /* initial */

      if ( fnd == 1 )             /* interpolate obs_tgain using evt_mjd_obs */
      {
         t_yi_ww1 = inp_p->tgainVal[tgIdx(yi,ww+1)];     /* tgain[yi,ww1] */
         
        /* dph spec eq.3  */
         if ( hpe_gt(timeGrid[1], timeGrid[0])==1 )         /* ie. a2 >0 */
         {
            inp_p->obs_tgain[yi] = t_yi_ww0 + (a1/a2)*(t_yi_ww1 - t_yi_ww0); 
         }
      } /* end: if fnd==1 */
   }

  /* ----------------------------------------------------
   * dph spec step4:  compute gainmap_2nd_order (eg.4 )
   *             (ie. compute G_2nd[ xj=0:47, yi=0:575] )
   * ----------------------------------------------------*/

   int RAWX_LEN  =  ( int ) inp_p->rawxgridSize ;   /* 8/2012 */

   inp_p->G_2nd=(double*)calloc(RAWX_LEN*RAWY_LEN,sizeof(double));
   long xj ;
   for (xj=0; xj<RAWX_LEN; xj++)        /* 0->47 */
   {
      for (yi=0; yi<RAWY_LEN; yi++)     /* 0->575 */
      {
         {
            /* dph spec eq. (4)  */
            inp_p->G_2nd[g2Idx(xj,yi,RAWX_LEN)] =      /*8/2012-pass RAWX_LEN*/
                 inp_p->gainmapVal[gmIdx(1,xj,yi,RAWX_LEN)]/inp_p->obs_tgain[yi] ; 
         }
      }
   }

  /* ------------------
   * close the table 
   * ------------------*/
   dmTableClose( srcBlock );

   return ;
}  /* end: load_S_new_gain_table()   */

/*------------------------------------------------------------------------
 * check new gainFile's pi_double limits and set the status bit if pi>max 
 *
 * This function can be used for both new hrcS and new hrcI gain.
 *------------------------------------------------------------------------ */
void check_spi_limit(
               EVENT_REC_P_T evt_p /*U*/
                  )
{
   /* for new hrcS and new hrcI gain file, use pi_double to check the range,
    * then convert it to PI; 
    */
   if (evt_p->pi_double > HDET_MAX_PI_NEW )     /* 1023 */
   {
      evt_p->status |= HDET_PI_VALUE_STS;       /* set status bit7 for PI*/
      evt_p->pi = HDET_MAX_PI_NEW ;             /* 1023 */
   }
   else if ((evt_p->pi_double < DS_HRC_PI_TLMIN) ||
            ds_dNAN(evt_p->pi_double) )    /* PI=NAN should not occur*/
   {
      evt_p->pi = DS_HRC_PI_TLMIN;        /* 0 */
   }
   else
   {
      /* 11/2009 : change outCol PI dtype from long to short for new gain map.
       *           Here, we've evt_p->pi set as LONG dtype, 
       *           so, let's still keep pix_double_to_long().
       *   Then, write_hrc_events() will call dmSetScalar_s(.,(short)evt_p->pi);
       */
      /*evt_p->pi=(short)evt_p->pi_double; Linux: 11.8->11  "-13.7"->"-13" */
      pix_double_to_long(evt_p->pi_double,FALSE,&evt_p->pi);/*pi: 11.8->12; 0.7->1*/
   }

   return ;
} /* end: check_spi_limit() */

/*-----------------------------------------------------------------------------
 * calculate gain index for old 2dim gain image or for new hrcI 2dim gain image
 * (ie. gainflag=OLD_SI_GAIN || NEW_I_GAIN )  
 *
 * 10/2009-rename old_gain_index() to image_2dim_gain_index(); no code changes;
 *----------------------------------------------------------------------------- */
void image_2dim_gain_index(           
            INPUT_PARMS_P_T inp_p,    /*I: */
                EVENT_REC_T *evt_p    /*U: */
                   )
{
   /* JCC(7/18/00) - add 1 to gain_index */
   if (inp_p->gain_cdelt[0])
      evt_p->gain_index[0] = evt_p->rawpos[0] / inp_p->gain_cdelt[0] + 1 ;
   if (inp_p->gain_cdelt[1])
      evt_p->gain_index[1] = evt_p->rawpos[1] / inp_p->gain_cdelt[1] + 1 ;
} /* end: image_2dim_gain_index() */
