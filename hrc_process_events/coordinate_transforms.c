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
* FILE NAME: coordinate_transforms.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: 

  This file contains the functions used by hrc_process_events to transform  
  coarse position into various coordinate systems.


* NOTES:

  wmclaugh@cfa	Mar 28, 1996  First Version.

* REVISION HISTORY:
  JCC(5/1/00)-no code changes for the tap corrections, but add comments 
              in calc_coarse_coords().
  JCC(5/11/00)- use amps_dd for computation
  JCC(7/18/00)- add 1 to gain_index.
* (3/2005)-before l1h_coarse_to_chip, copy corr amp_sf from 
           evt to d_p and will be used for #3 degap.
* (2/2005)(6/2005)-remove unused variable
* (1/2009) - add hrcS 3dim gain image  ( obsolete 10/2009 )
* 10/2009 - replace peter hrcS 3dim gain image w/ dph new hrcS gain table.
* 10/2009 - add fap hrcI new gain image.
*H***********************************************************************/

#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

#ifndef L1_ASPECT_DEFS_H
#include "l1_aspect_defs.h"
#define L1_ASPECT_DEFS_H
#endif 



void calc_coarse_coords(
   EVENT_REC_T     *evt_p,  /* I/O structure containing event data       */
   INPUT_PARMS_P_T   inp_p, /* I  ptr to  struct containing input parms  */
   STATISTICS_T    *stat_p, /* O   structure containing statistical cnts */
   DEGAP_CONFIG_P_T d_p,    /* I/O - degap configuration structure       */
   dsErrList*       err_p)  /* I/O - error list                          */

{
   double    fine[HDET_NUM_PLANES];    /* initial fine position          */
   int       coarse[HDET_NUM_PLANES];  /* coarse position                */
   int       plane;
   boolean   bad_denom = FALSE; 

   /* get pha sum */
   sum_phas_hrc(evt_p);  /*JCC: return 'evt_p->amp_tot[plane] and sum_amps' */

   if (inp_p->gainflag == NEW_S_GAIN )   /* 10/2009 */
   {
      calc_DDn( evt_p ) ;           /* compute DDn */
   }

   for (plane = HDET_PLANE_X; plane < HDET_NUM_PLANES; plane++)
   {
      /* check for wire charge in plane */
      if ((inp_p->wire_charge == HDET_WIRE_ON) &&
         ((evt_p->amps_dd[plane][HDET_3RD_AMP] >
           evt_p->amps_dd[plane][HDET_2ND_AMP]) ||
           (evt_p->amps_dd[plane][HDET_1ST_AMP] > 
           evt_p->amps_dd[plane][HDET_2ND_AMP])))
      {
         stat_p->bad_dist[plane]++;

         if (plane == HDET_PLANE_X)
         {
            evt_p->status |= HDET_U_CNTR_STS; 
         }
         else
         {
            evt_p->status |= HDET_V_CNTR_STS; 
         } 
      } 

      /*-------------------------------------------------------------
       * JCC(5/1/00) - notes :
       *    'check_tap_ring' returns a NEW amps_dd[][3RD_AMP]
       *    'sum_phas_hrc()' returns (amps_dd[1] + amps_dd[2]+ tap_ring_A3
       *------------------------------------------------------------*/
      if (evt_p->amp_tot[plane] > 0)
      {

         fine[plane] = (evt_p->amps_dd[plane][HDET_3RD_AMP] -
            evt_p->amps_dd[plane][HDET_1ST_AMP]) / evt_p->amp_tot[plane];

      }
      else
      {
         fine[plane] = 0; 
         bad_denom = TRUE; 
         evt_p->status |= HDET_ZERO_SUM_STS;
      }
      coarse[plane] = evt_p->cp[plane]; 
   } /* end for plane */  


   l1h_coarse_to_raw(coarse, fine, evt_p->rawpos);/*recompute evt_p->rawpos[x,y]*/

   d_p->amp_sf = evt_p->amp_sf ;   /* for degap #3 */
   l1h_coarse_to_chip(d_p, coarse, fine, evt_p->chippos, &evt_p->chipid,err_p);  

   if (inp_p->gainflag == NEW_S_GAIN)  
   {
     /* ---------------------------------------------------------------------- 
      * 10/2009 - For new hrcS gain table, use evt_p->rawpos to get gain index
      * and pi_double. PI value will be readjusted in calculate_pi_hrc().
      * ----------------------------------------------------------------------*/
      S_new_gain_index_pi(inp_p, evt_p); /*load gainTab data + compute pi_double*/
   }
   else   /* 10/2009 - for both OLD_SI_GAIN && NEW_I_GAIN */
   {
      image_2dim_gain_index(inp_p, evt_p);  /* load gain image data (only!)*/
     /* pi||pi_double will be computed in calculate_pi_hrc() later.*/
   }

   /* update dependency masks to include newly calculated fields */
   stat_p->dependencies |= HDET_MASK_CHIP_REQ;  

   /* perform ratio checks */
   if (inp_p->do_ratio)
   {
      ratio_checks_hrc(evt_p, inp_p, stat_p);
   }

   if (bad_denom)
   {
      stat_p->bad_bot++;
      evt_p->status |= HDET_FIN_POS_STS; 
   }
} /* end: calc_coarse_coords() */



boolean calc_chip_coords(
   EVENT_REC_P_T   evt_p,  /* I/O structure containing event data        */
   INPUT_PARMS_P_T inp_p,  /* I  ptr to  struct containing input parms   */
   STATISTICS_P_T  stat_p, /* O   structure containing statistical counts*/
   DEGAP_CONFIG_P_T d_p,   /* I/O - degap configuration structure        */
   dsErrList*       err_p) /* I/O - error list                           */
{
   boolean err = FALSE; 

   switch(inp_p->start)
   {
      case HDET_COARSE_VAL:
         calc_coarse_coords(evt_p, inp_p, stat_p, d_p, err_p);
      break; 

      case HDET_CHIP_VAL:
         if (inp_p->hrc_system == HRC_IMAGE_INST)
         {
            evt_p->chipid = 0;
         }
      break; 

      case HDET_TDET_VAL:
         pix_tdet_to_chip(evt_p->tdetpos, &evt_p->chipid, evt_p->chippos);

         /* temporary fix */ 
         if (inp_p->hrc_system == HRC_IMAGE_INST)
         {
            evt_p->chipid = 0;
         }

      break; 

      case HDET_DET_VAL: /* fallthrough intended */ 
      case HDET_TAN_VAL: /* fallthrough intended */ 
      case HDET_SKY_VAL: /* fallthrough intended */ 
      case HDET_NONE_VAL: /* fallthrough intended */  
      default: 
         /* invalid */ 
         err = TRUE; 
      break; 

   }
  
   return (err); 
}




boolean calculate_coords_hrc(
   EVENT_REC_P_T   evt_p,   /* I/O structure containing event data         */
   INPUT_PARMS_P_T inp_p,   /* I  ptr to  struct containing input parms    */
   STATISTICS_P_T  stat_p,  /* O   structure containing statistical counts */
   ASPECT_REC_P_T  aspect,  /* I   aspect information                      */
   DEGAP_CONFIG_P_T d_p,    /* I/O - degap configuration structure         */
   short asp_type_flag,     /* I   what type of aspect correction is it?   */
   dsErrList*       err_p)  /* I/O - error list                            */
{
   boolean err = FALSE;  /* false = no error occurred */  

#ifdef NOT_MOVED_FOR_ADC_CORR 
   /* if HRC-i flight data extra bit should be removed */ 
   if ((evt_p->cp[HDET_PLANE_Y] >= 64) && 
       (inp_p->hrc_system == HRC_IMG_SYS)) 
   {
      evt_p->cp[HDET_PLANE_Y] -= 64; 
   } 
#endif
  
   if (((evt_p->cp[HDET_PLANE_X] < d_p->min_tap[HDET_PLANE_X]) || 
       (evt_p->cp[HDET_PLANE_X] > d_p->max_tap[HDET_PLANE_X])) || 
       ((evt_p->cp[HDET_PLANE_Y] < d_p->min_tap[HDET_PLANE_Y]) || 
       (evt_p->cp[HDET_PLANE_Y] > d_p->max_tap[HDET_PLANE_Y])))   
   { 
      err = TRUE;
   }
   else
   { 
      if (!(calc_chip_coords(evt_p, inp_p, stat_p, d_p, err_p)))
      {
         /* register long chipint;  */
         evt_p->workpos[HDET_PLANE_X] = evt_p->chippos[HDET_PLANE_X]; 
         evt_p->workpos[HDET_PLANE_Y] = evt_p->chippos[HDET_PLANE_Y]; 

         if (inp_p->randpixsize > FLT_EPSILON)
         {
            double rand_x, rand_y;
   
            /* apply a random value from -randpixsize to +randpixsize
             * to avoid quantization artifacts. (pix_short_to_double
             * doesn't accept 0 is input so use 1 and subtract it from
             * the returned value.
             */
            pix_short_to_double(1, 1, &rand_x);
            pix_short_to_double(1, 1, &rand_y);
            rand_x -= 1.0;
            rand_y -= 1.0;
            rand_x *= (2.0 * inp_p->randpixsize);
            rand_y *= (2.0 * inp_p->randpixsize);
            evt_p->workpos[HDET_PLANE_X] += rand_x; 
            evt_p->workpos[HDET_PLANE_Y] += rand_y; 
         }

         switch(inp_p->stop)
         {
            case HDET_SKY_VAL:  /* fallthrough intended */  
               pix_chip_to_fpc(evt_p->chipid, evt_p->workpos, evt_p->fppos);
 
               if (inp_p->processing == HRC_PROC_FLIGHT)
               {
                  VEC2_DBLE fpc;
   
                  fpc[0] = evt_p->fppos[HDET_PLANE_X]; 
                  fpc[1] = evt_p->fppos[HDET_PLANE_Y]; 

		  if (asp_type_flag == ASP_FTYPE_OFFSETS)
		  {  
		     pix_apply_aspect(fpc, aspect->asp_sol, evt_p->skypos); 
		  }
		  else
		  {
		     double cel[2];
		     
		     dmTanPixToWorld(fpc, aspect->asp_sol, inp_p->crpix,
                               inp_p->cdelt, cel);
		     dmTanWorldToPix(cel, inp_p->crval, inp_p->crpix,
                               inp_p->cdelt, evt_p->skypos);
   		  }
               }
               else  
               {
                  evt_p->skypos[HDET_PLANE_X] = evt_p->fppos[HDET_PLANE_X];
                  evt_p->skypos[HDET_PLANE_Y] = evt_p->fppos[HDET_PLANE_Y];
               }

            case HDET_DET_VAL:  /* fallthrough intended */ 
               if (inp_p->processing == HRC_PROC_FLIGHT)
               {
                  if (inp_p->stop != HDET_SKY_VAL)
                  {
                     pix_chip_to_fpc(evt_p->chipid, evt_p->workpos, 
                                     evt_p->fppos);
                  } 
                  evt_p->detpos[HDET_PLANE_X] = evt_p->fppos[HDET_PLANE_X];
                  evt_p->detpos[HDET_PLANE_Y] = evt_p->fppos[HDET_PLANE_Y];
               }
               else
               {
                  pix_chip_to_fpc(evt_p->chipid, evt_p->workpos, 
                                  evt_p->detpos);
               } 
            case HDET_TDET_VAL: /* fallthrough intended */ 
               pix_chip_to_tdet(evt_p->chipid, evt_p->chippos,
                                evt_p->tdetpos);
            break; 
            
            default:
            break; 
         }
      }
   }

   return (err); 
}

short parse_coord_range(char* string) /* I - start/stop string */
{
   short val = HDET_INV_VAL;
 
   if (strcmp(string, HDET_NO_COORD) == 0)
   {
      val = HDET_NONE_VAL;
   }
   if (strcmp(string, HDET_COARSE_COORD) == 0)
   {
      val = HDET_COARSE_VAL;
   }
   else if (strcmp(string, HDET_CHIP_COORD) == 0)
   {
      val = HDET_CHIP_VAL;
   }
   else if (strcmp(string, HDET_TDET_COORD) == 0)
   {
      val = HDET_TDET_VAL;
   }
   else if (strcmp(string, HDET_DET_COORD) == 0)
   {
      val = HDET_DET_VAL;
   }
   else if (strcmp(string, HDET_TAN_COORD) == 0)
   {
      val = HDET_TAN_VAL;
   }
   else if (strcmp(string, HDET_SKY_COORD) == 0)
   {
      val = HDET_SKY_VAL;
   }
 
   return (val);
}


 
boolean map_start_column (
   short  start_pos,      /* I  - starting position of transformations */
   short* x_col,          /* O  - column mapping for x axis */
   short* y_col)          /* O  - column mapping for y axis */
{
   boolean err = FALSE;
 
   switch(start_pos)
   {
      case HDET_NONE_VAL:     /* fallthrough intended */
      case HDET_COARSE_VAL: /* fallthrough intended */ 
      case HDET_CHIP_VAL:
         *x_col = HDET_CHIP_X;
         *y_col = HDET_CHIP_Y;
      break;
      case HDET_TDET_VAL:
         *x_col = HDET_TDET_X;
         *y_col = HDET_TDET_Y;
      break;
      case HDET_DET_VAL:
         *x_col = HDET_DET_X;
         *y_col = HDET_DET_Y;
         err = TRUE;
      break;
      case HDET_TAN_VAL:
         *x_col = HDET_TAN_X;
         *y_col = HDET_TAN_Y;
         err = TRUE;
      break;
      case HDET_SKY_VAL:
         *x_col = HDET_SKY_X;
         *y_col = HDET_SKY_Y;
         err = TRUE;
      break;
      default:
         err = TRUE;
      break;
   }
 
   return (err);
}
