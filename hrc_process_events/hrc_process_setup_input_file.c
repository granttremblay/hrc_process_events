/*                                                                
**  Copyright (C) 1998-2009  Smithsonian Astrophysical Observatory 
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
 
* FILE NAME: hrc_process_setup_input_file.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: The routine hrc_process_setup_input_file is responsible for
  opening the input event file. It takes in a pointer to a structure of
  INPUT_PARMS_T which contains all of the relevant parameter data needed 
  from the input parameter file to open the event file. The input file
  housekeeping information is then written to the EVENT_SETUP_T struct
  passed into the function. The routine does not return a value. If any 
  errors are detected they are stored in an error list and retrieved 
  latter via the error lib by the calling routine.
 
JCC(8/2002) - get obs_tstop from time_stop (time_stop is a string of 'TSTOP').
JCC(11/2002)- if ra/dec_nom not get from obsfile, then get them from evt1 infile;
              else ( make sure obsfile and evt file have same value);
(2/2005)-add dmBlockGetNo; remove dmDatasetGetKernel;
(1/2009)-infile-stack and infile-access will be checked in hpe.c 
        -add hpePrintErr()
10/2009-get mjd_obs from evtfile for dph new hrcS gain table. 
       -rename MAX_INST_KYWRD_LEN to HPE_LEN_80.
*H***********************************************************************/

#include <stdio.h>
#include <ctype.h>
 
#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#define HRC_PROCESS_EVENTS_H
#endif

void hrc_process_setup_input_file(
   EVENT_SETUP_P_T evtin_p,   /* I/O - input event file housekeeping ptr*/
   INPUT_PARMS_P_T inp_p,     /* I/O - input parameter/run data pointer */ 
   STATISTICS_P_T  stat_p,    /* O   - statistics data structure ptr    */
   dsErrList*      hpe_err_p) /* O   - error list pointer               */
{
   short   rr, cc; 
   char**  in_names;       /* input evt file col names    */
   boolean alloc_failure = FALSE; 

   /* update statistical file counts */
   stat_p->num_files_in++;

  /* ------------------------------------------------------------------
   * 1/2009: assume "infile-stack & infile-access were already checked 
   *         in hpe.c.  So, we won't need to verify them again here.  
   * ------------------------------------------------------------------ */

      long           temp_cols;  /* count of columns (vectors count as 1) */
      long*          temp_dim  = NULL; /* array of column dimensions      */
      dmDescriptor** temp_desc = NULL; /* array of column/vector desc's   */

      ds_open_event_file(evtin_p->file, &evtin_p->dataset,
                         &evtin_p->extension);
 
      if (dmBlockGetNo(evtin_p->extension) > 1 )
         evtin_p->primary = dmDatasetMoveToBlock(evtin_p->dataset, 1);
      else
         evtin_p->primary = NULL ;
 
      temp_cols = dmTableGetNoCols(evtin_p->extension);
      evtin_p->num_cols = 0;
      evtin_p->dim = NULL; 
 
      if (((temp_dim = (long*) calloc(temp_cols, sizeof(long))) != NULL) &&
          ((temp_desc = (dmDescriptor**)
          calloc(temp_cols, sizeof(dmDescriptor*))) != NULL))
      {
         for (rr = 0; rr < temp_cols; rr++)
         {
            /* get column or vector descriptors */
            temp_desc[rr] = dmTableOpenColumnNo(evtin_p->extension, rr+1);

            /* get the total number of column components */ 
            evtin_p->num_cols +=
               (temp_dim[rr] = dmGetElementDim(temp_desc[rr]));
         }
      }
      else
      {
         alloc_failure = TRUE;
      }

      /* allocate necessary memory */
      evtin_p->mapping = (short*) calloc(evtin_p->num_cols, sizeof(short));
      evtin_p->desc = (dmDescriptor**) calloc(evtin_p->num_cols,
         sizeof(dmDescriptor*));
      evtin_p->types = (dmDataType*) calloc(evtin_p->num_cols,
         sizeof(dmDataType));
      if ((in_names = (char**) calloc(evtin_p->num_cols,
         sizeof(char*))) != NULL)
      {
         for (rr = 0; rr < evtin_p->num_cols; rr++)
         {
            if ((in_names[rr] = (char*) calloc(DS_SZ_PATHNAME, sizeof(char))) ==
                 NULL)
            {
               alloc_failure = TRUE;
            } 
         }  
      }  
 
      if ((evtin_p->desc == NULL) || (in_names == NULL) ||
         (evtin_p->mapping == NULL) || alloc_failure)
      {
         dsErrAdd(hpe_err_p, dsALLOCERR, Individual, Custom,
            "ERROR: Memory allocation failed setting up input event file.");
      }
      else
      {
         char    telescop[HPE_LEN_80];
         char    datamode[HPE_LEN_80];
         char    detnam[HPE_LEN_80]; 
         double  time_val;
         double  evt_ra_nom, evt_dec_nom ;   /* 11/2002 */

         for (cc = rr = 0; rr < temp_cols; rr++)
         {
            if (temp_dim[rr] == 1)
            {
               /* column is a non-vector column - get info */ 
/****
               evtin_p->desc[cc] =
                  dmTableOpenColumnNo(evtin_p->extension, rr+1);
****/
               memcpy(&evtin_p->desc[cc], &temp_desc[rr], 
                      sizeof(dmDescriptor*));

               evtin_p->types[cc] = dmGetDataType(evtin_p->desc[cc]);
               dmGetName(evtin_p->desc[cc], in_names[cc], DS_SZ_PATHNAME);
               cc++;
            }
            else
            {
               long zz;
 
               /* column is a vector column - get info for each component */ 
               for (zz = 0; zz < temp_dim[rr]; zz++)
               {
                  evtin_p->desc[cc] =
                     dmGetCpt(temp_desc[rr], zz+1);
                  evtin_p->types[cc] = dmGetDataType(evtin_p->desc[cc]);
                  dmGetName(evtin_p->desc[cc], in_names[cc], DS_SZ_PATHNAME);
                  cc++;
               }
            }
         }
 
         /* free memory used to obtain vector components */
         if (temp_desc != NULL)
         {
            free(temp_desc);
            temp_desc = NULL; 
         }
         if (temp_dim != NULL)
         {
            free(temp_dim);
            temp_dim = NULL; 
         }

         /* set up mappings */
         parse_hrc_evt_columns(in_names, evtin_p->num_cols, evtin_p->mapping); 

         for (rr = evtin_p->num_cols; rr--; )
         {
            inp_p->scl_xsts |= (evtin_p->mapping[rr] == HDET_AMP_SF);
            if (in_names[rr])
            {
               free(in_names[rr]);
            } 
         }

/****
         if (dmKeyRead_d(evtin_p->extension, inp_p->time_start, &time_val) 
             != NULL)
         {
            if (time_val >= inp_p->end_time)
            {
               inp_p->start_time = time_val;
            } 
            else 
            {
               dsErrAdd(hpe_err_p, dsHPETIMEORDERERR, Individual, Generic);
            }
         }
         inp_p->end_time = inp_p->start_time;
         if (dmKeyRead_d(evtin_p->extension, inp_p->time_stop, &time_val) 
             != NULL)
         {
            if (time_val >= inp_p->end_time)
            {
               inp_p->end_time = time_val;
            }
            else 
            {
               dsErrAdd(hpe_err_p, dsHPETIMEORDERERR, Individual, Generic);
            }
         }
****/
         if ((inp_p->use_obs == 0) &&
             (dmKeyRead_d(evtin_p->extension, inp_p->time_start, &time_val)
              != NULL))
         {
            inp_p->obs_tstart = time_val; 
         } 
         if ((inp_p->use_obs == 0) &&
             (dmKeyRead_d(evtin_p->extension, inp_p->time_stop, &time_val)
              != NULL))
         {
            inp_p->obs_tstop = time_val; 
         } 

         if ((inp_p->use_obs == 0) && 
             (dmKeyRead_c(evtin_p->extension, TELESCOP_KEY, telescop, 
             HPE_LEN_80) != NULL))
         {
            short ii; 

            if ((ii = strlen(telescop)) > HPE_LEN_80)
            {
               ii = HPE_LEN_80 - 1;
            }
            telescop[ii] = '\0';
            while (ii-- > 0)
            {
               telescop[ii] = (char) toupper(telescop[ii]);
            }

            if ((strcmp(telescop, TELESCOP_XRCF_KEY) == 0) ||
                (strcmp(telescop, TELESCOP_XRCF_KEY2) == 0))
            {
               inp_p->processing = HRC_PROC_XRCF;
            }
         } 

         /* get DATAMODE keyword to see if next-in-line mode */
         if ((inp_p->use_obs == 0) && 
             (dmKeyRead_c(evtin_p->extension, DATAMODE_KEY, datamode,
              HPE_LEN_80) != NULL))
         {
            inp_p->next_in_line = 
               ((ds_strcmp_cis(datamode, NEXT_IN_LINE_VAL) == 0) ||  
                (ds_strcmp_cis(datamode, NEXT_IN_LINE_VAL2) == 0)); 
         }

         /* get DETNAM from input file header if no obs.par read */ 
         if ((inp_p->use_obs == 0) && 
             (dmKeyRead_c(evtin_p->extension, "DETNAM", detnam,
              HPE_LEN_80) != NULL))
         {
            strcpy(inp_p->instrume, detnam); 
         }

        /*------------------------------------------
         *(11/2002) get RA_NOM from evt1 hdr 
         *-----------------------------------------*/
         if (dmKeyRead_d(evtin_p->extension, "RA_NOM", &evt_ra_nom) != NULL)
         {
            if (inp_p->ra_nom == NOM_NOT_FOUND)
            {
               inp_p->ra_nom =  evt_ra_nom ;
            }
            else
            {
                 if ( inp_p->ra_nom != evt_ra_nom )
                 {
                    dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                    "WARNING: ra_nom in obsfile and infile not same. taking from obsfile.\n");
                 }
            } /* if (inp_p->ra_nom == NOM_NOT_FOUND ) */
         }
         else
         {
            /* not found in infile and obsfile */
            if (inp_p->ra_nom==NOM_NOT_FOUND)
            {
               dsErrAdd(hpe_err_p, dsFINDKEYWORDFERR, Individual, Custom,
                  "ERROR : ra_nom is required in obsfile or infile.\n");
            }
         } /* if dmKeyRead_d(evtin,RA_NOM) */

        /*------------------------------------------
         *(11/2002) get DEC_NOM from evt1 hdr 
         *-----------------------------------------*/
         if (dmKeyRead_d(evtin_p->extension, "DEC_NOM", &evt_dec_nom) != NULL)
         {
            if (inp_p->dec_nom == NOM_NOT_FOUND)
            {
               inp_p->dec_nom =  evt_dec_nom ;
            }
            else
            {
                if ( inp_p->dec_nom != evt_dec_nom )
                {
                   dsErrAdd(hpe_err_p, dsGENERICERR, Individual, Custom,
                  "WARNING: dec_nom in obsfile and infile not same. taking from obsfile.\n");
                }
            } /* if (inp_p->dec_nom == NOM_NOT_FOUND ) */
         }
         else
         {
            /* not found in infile and obsfile */
            if (inp_p->dec_nom==NOM_NOT_FOUND)
            {
               dsErrAdd(hpe_err_p, dsFINDKEYWORDFERR, Individual, Custom,
                    "ERROR : dec_nom is required in obsfile or infile.\n");
            }
         } /* if dmKeyRead_d(evtin,DEC_NOM) */

        /*----------------------------------------------------------------
         *(10/2009) get event key MJD_OBS for new hrcS gain table 
         * If it is missing and has no enough info to compute, Error out!
         *----------------------------------------------------------------*/
         inp_p->mjd_obs_warn = 0;
         if (dmKeyRead_d(evtin_p->extension, "MJD_OBS", &inp_p->evt_mjd_obs) == NULL)
         {
            inp_p->evt_mjd_obs = 0.0 ; 

           /* issue a warning later for new hrcS gain table */
            inp_p->mjd_obs_warn = 1; /*ERROR OUT: mjd_obs isn't in evtfile & can't be computed*/


           /* compute evt_mjd_obs :  MJD_OBS=MJDREF+(TIMEZERO+TSTART)/86400.0 */
            double tmp_mjdref = 0.0 ;
            double tmp_timezero= 0.0 ;
            double tmp_tstart= 0.0 ;
            if ((dmKeyRead_d(evtin_p->extension,"MJDREF",&tmp_mjdref)==NULL) ||
                (dmKeyRead_d(evtin_p->extension,"TIMEZERO",&tmp_timezero)==NULL)||
                (dmKeyRead_d(evtin_p->extension,"TSTART",&tmp_tstart)==NULL)
               )
            {
               /*ERROR OUT: mjd_obs isn't in evtfile and no enough info to compute it */
               /*here let's still keep  inp_p->evt_mjd_obs=0.0 and inp_p->mjd_obs_warn=1*/
            }
            else
            {
              /*mjd_obs isn't in evtfile but can be computed */
               inp_p->evt_mjd_obs = tmp_mjdref + (tmp_timezero+tmp_tstart)/86400.0;
               inp_p->mjd_obs_warn = 2; /*WARNING: mjd_obs is computed; no ERROR OUT*/
            }
         }

         /* set data dependency mask */ 
         stat_p->dependencies = dependency_check_init(evtin_p->mapping,
                                                      evtin_p->num_cols);

         /* determine if data dependencies have been met */
         dependency_check_hrc(inp_p, stat_p, hpe_err_p);

         free(in_names); 
      }

} /* end: hrc_process_setup_input_file() */


/* -----------------------------------------------
 *   1/2009 : new routine to display error messge 
 * ----------------------------------------------- */
dsErrCode hpePrintErr( 
                dsErrList*  hpe_err_p, /* U */
                  FILE*     log_ptr,   /* I */
                   int      debug      /* I */
                     )
{
   dsErrCode err = dsNOERR;

   if (hpe_err_p->size > 0)
   {
      if (hpe_err_p->contains_fatal != 0)
         err = dsGENERICERR;

      dsErrPrintList(hpe_err_p, dsErrTrue);

      if ((debug > DEBUG_LEVEL_0) && (log_ptr != stderr) &&
            (log_ptr != stdout))
      {
         dsErrDirectOutput(log_ptr);
         dsErrPrintList(hpe_err_p, dsErrTrue);
         dsErrDirectOutput(stderr);
      }

      while (hpe_err_p->size > 0) dsErrRemove(hpe_err_p);
   }

   /*printf("Fatal=%ld, errNum=%ld\n",hpe_err_p->contains_fatal,hpe_err_p->size);*/

   return err ;

} /* end: hpePrintErr() */
