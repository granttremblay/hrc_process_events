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

/*H**************************************************************************
 
* FILE NAME: calculate_pi_hrc.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: 

  This routine will calculate the pulse invarience for a given event. The
  value will be passed back via a pointer to a structure containing all
  relevant event data (EVENT_REC_T). The routine is called by hrc_process_events.

  The computational algorithm for pulse invarience has yet to be determined. 
 
* NOTES:

  wmclaugh@cfa	Mar 28, 1996  First Version.

* JCC(6/25/01) - if (PI > 255), set to 255.  (for old gain image)

1/2009 - add peterR hrcS 3dim gain image.    ( 10/2009: obsolete )
10/2009- replace peterR hrcS 3dim gain image w/ dph new hrcS gain table.
       - 5 columns in the new hrcS gain table : 
               1   GAINMAP[2,48,576]  real8
               2   TGAIN[576,18]      real8
               3   RAWXGRID[48]       real8
               4   RAWYGRID[576]      real8
               5   TIMEGRID[18]       real8
10/2009- add fap new hrcI gain image which has the key sampnorm.
       - see 'Notes on outCol PI'
*H**************************************************************************/

#include "hrc_process_events.h"
#include "dslib.h"

void calculate_pi_hrc(
   float*        gain_p,   /* I   gain image data                       */
   INPUT_PARMS_P_T inp_p,  /* I   input parms (gain image axes lengths) */
   EVENT_REC_P_T evt_p)    /* I/O structure containing event data       */ 
{

if (inp_p->gainflag == OLD_SI_GAIN )    /* old 2dim gain image*/
{
   long offset = (evt_p->gain_index[1] -1) * inp_p->gain_axlen[0] + 
                 evt_p->gain_index[0] -1;   

   /* PHA value * (factor from gain image) */ 
   evt_p->pi = evt_p->pha * gain_p[offset]; 

   /*-- old gain img:  use PI to check the range --- */
   if (evt_p->pi > HDET_MAX_PI_OLD )     /* 255 */
   {
      evt_p->status |= HDET_PI_VALUE_STS;
      evt_p->pi = HDET_MAX_PI_OLD ;     /* 255 */
   }
}
else if (inp_p->gainflag == NEW_I_GAIN )    /* 10/2009-new hrcI gain image*/
{
   long offset = (evt_p->gain_index[1] -1) * inp_p->gain_axlen[0] +
                 evt_p->gain_index[0] -1;

   /* fap: samp=evt_p->sumamps*2**(evt_p->amp_sf-1)/C148_from_sampnorm_key;
           evt_p->pi = samp * gain_p[offset];    */
    double samp = ldexp(evt_p->sum_amps,evt_p->amp_sf-1)/inp_p->sampnorm; 
    evt_p->pi_double = samp * gain_p[offset] ;

  /* 10/2009 - set pi_double min/max properly */
    check_spi_limit( evt_p ) ;
}
else if (inp_p->gainflag == NEW_S_GAIN)   /*10/2009-new hrcS gain table*/
{
  /*(10/2009)NOTE: pi_double was computed earlier in S_new_gain_index_pi()*/

  /* 10/2009 - set pi_double min/max properly */
   check_spi_limit( evt_p ) ;
}
else
{
  /* shouldn't happen */
}

} /* end: calculate_pi_hrc() */

/* ---------------------------------------------------------------------- 
 * load gain data for old 2dim gain image or new hrcI 2dim gain image or 
 * new hrcS gain table.
 * ----------------------------------------------------------------------*/

void load_gain_image(
   char*           filename,     /* I   - name of gain file          */
   INPUT_PARMS_P_T inp_p,        /* O   - input parms (store cdelts) */
   float**         gain_p,       /* O   - pointer to gain map        */ 
   dsErrList*      err_p)        /* O   - error list                 */
{
   /* initialize */
   inp_p->gainflag = OLD_SI_GAIN; 
   inp_p->gain_axlen[0] = inp_p->gain_axlen[1] = 0;

   if ((ds_strcmp_cis(filename, "NONE") == 0) || 
       (ds_strcmp_cis(filename, "\0") == 0))
   {
      /* no gain map- do nothing. no error out; */ 
      /* printf("gainfile is NONE or blank(%s)\n",inp_p->gain_file);*/
   } 
   /* else if (dmDatasetAccess(filename,"R")==dmTRUE) */
   else if (dmFileExists(filename)==DMF_NOFILE) 
   {
      /* file is not NONE but can't find it - error out*/ 
      dsErrAdd (err_p, dsOPENFILEFERR, Individual, Custom, "ERROR: can't find gain file, %s\n", inp_p->gain_file );  /*not dsGENERICERR*/
      return ;
   }
   else    /* find gain file */ 
   {
     /* ---------------------------------------------
      * set gainflag to NEW_S_GAIN or OLD_SI_GAIN.
      * --------------------------------------------*/
      set_gainflag (inp_p, err_p ) ;       
      /* printf("inp_p->gainflag = %d\n", inp_p->gainflag);*/

     /* ----------------------------------------------------
      *    BEGIN: load data for new hrcS gain table & STOP !
      * ---------------------------------------------------- */
      if ( inp_p->gainflag == NEW_S_GAIN )    
      {
         load_S_new_gain_table( inp_p ) ; /* load gain data */
         return ;                         /* STOP! */
      }
     /* --------------------------------------------------
      *    END: load data for new hrcS gain table & STOP ! 
      * -------------------------------------------------- */

     /* -----------------------------------------------------------------
      * BEGIN - load data for old 2dim gainImg OR new hrcI 2dim gainImg.
      *
      * NOTE - old 2dim gain image and new hrcI 2dim gain image have same
      *        data format except SAMPNORM only exists in new hrcI gain.
      * -----------------------------------------------------------------*/
      long    image_area;          /* size of image (length*width) */ 
      long*   axes_data;           /* array for image dimensions   */ 
      dmBlock* file_desc;          /* file block pointer           */ 
      dmDescriptor* data_desc;     /* image dsescriptor            */  

      file_desc = dmImageOpen(filename);
      data_desc = dmImageGetDataDescriptor(file_desc);

     /* -----------------------------------------------------------------
      * Here we starts w/  gainflag=OLD_SI_GAIN. Further check on the
      * SAMPNORM key in the gain image hdr. If it's found, then
      * change gainflag to NEW_I_GAIN (ie. it's fap new hrcI gain image).
      * -----------------------------------------------------------------*/
      if ((dmKeyRead_d(file_desc, "SAMPNORM", &inp_p->sampnorm) != NULL))
      {
            inp_p->gainflag = NEW_I_GAIN ; 
      }
     /*-- end--*/

      /* axes_data[2] required for old gain image and new hrcI gain image */
      dmGetArrayDimensions(data_desc, &axes_data);

      double  crval[2],
              cdelt[2] = {0, 0};   /* wcs values for image         */
      float   crpix[2];            /* wcs values for image         */
      dmDescriptor* ax_grp[2];     /* image axis group descriptors */
      dmDescriptor* wcs_desc[2];      /* wcs coordinate descriptors   */

      /* compute image size */
      if ((image_area = (inp_p->gain_axlen[0] = axes_data[0]) * axes_data[1]) > 0)
      {
         if ((*gain_p = (float*)calloc (image_area,
            sizeof(float))) != NULL)
         {
            long img_LL[2] = {1,1}; /* lower left corner of image */
            long img_UR[2]; /* upper right corner of image */

            img_UR[0] = axes_data[0];
            img_UR[1] = axes_data[1];
       
            /* get data */
            dmImageDataGetSubArray_f(data_desc, img_LL, img_UR,  
                                     *gain_p);
         }
         else
         {
            /* error- allocation failed */ 
         }
      }
      else
      {
         /* error- bad image dimensions */ 
      }

      /* get CDELT information since its needed to compute gain indices */
      if ((ax_grp[0] = dmArrayGetAxisGroup(data_desc, 1)) != NULL) 
      {
         if ((wcs_desc[0] = dmDescriptorGetCoord(ax_grp[0])) != NULL)  
         {
            dmCoordGetTransform_f(wcs_desc[0], crpix, crval, cdelt, 1); 
            if (cdelt[0] > 0)  /* avoid division by zero and negatives */ 
            {
               inp_p->gain_cdelt[0] = cdelt[0]; 
            } 
         } 
      } 
      if ((ax_grp[1] = dmArrayGetAxisGroup(data_desc, 2)) != NULL) 
      {
         if ((wcs_desc[1] = dmDescriptorGetCoord(ax_grp[1])) != NULL) 
         {
            dmCoordGetTransform_f(wcs_desc[1], &crpix[1], &crval[1], 
                                  &cdelt[1], 1);
            if (cdelt[1] > 0) /* avoid division by zero and negatives */
            {
               inp_p->gain_cdelt[1] = cdelt[1];
            }
         }
      }

      free(axes_data);
      dmImageClose(file_desc);

     /* ----------------------------------------------------------------
      * END - load data for old 2dim gain image OR new 2dim gain image
      * ----------------------------------------------------------------*/
   } /*end:  find gain file */

} /* end: load_gain_image() */
