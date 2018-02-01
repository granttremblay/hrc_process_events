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

/*H***********************************************************************
 
* FILE NAME: adc_routines.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: 

  The file adc_routines.c contains the following modules used
  by hrc_process_events() in dealing with adc_correction tables:

        adc_table_load()  
        allocate_adc_table()  
        apply_adc_correction() 
        deallocate_adc_table()  
        load_adc_entry() 
        parse_adc_column() 
 
* NOTES:

  wmclaugh@cfa	July 14, 1998  First Version.

* REVISION HISTORY:
  JCC(5/11/00) - use amps_dd for computation
  JCC(2/9/01) - enhanced to memset axis[2] before calling dmGetScalar_c.
(2/2005)-initialize some variables.
*H***********************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

#ifndef ADC_COOR_DEFS_H
#include "adc_corr_defs.h"
#endif 


/* local function to parse column names from adc table */
static short parse_adc_column(char*);
 
/* local function to load entries from the adc table into dynamic memory */
static void load_adc_entry (ADC_SETUP_P_T, INPUT_PARMS_P_T,
                            ADC_CORR_P_T, ADC_CORR_P_T, dsErrList*);
 


/*H***********************************************************************

* DESCRIPTION: 

  The function adc_table_load is called by hrc_process_events to copy 
  the adc correction information from the user specified ADC file into
  memory. These values are later used to correct the raw adc values 
  in the input event file. 

  This routine opens the fits format adc file and reads in all of the 
  data from the ADCCORR extension of the file. It checks to see if the
  expected number of entries (corresponding to tap positions) have been
  read in. The routine expects that the calling routine has already 
  allocated sufficient memory to store the retrieved data. 
 
* NOTES:
 
The routine does not return n error status, but encountered errors are 
pushed onto the error list. 
 
*H***********************************************************************/

void adc_table_load(
    INPUT_PARMS_P_T  inp_p,     /* I   pointer to input parameters  */
    ADC_CORR_P_T     adc_x,   /* O   array of x coord adc corrections */ 
    ADC_CORR_P_T     adc_y,   /* O   array of y coord adc corrections */
    dsErrList*       err_p)   /* O   error list pointer               */
{
    short ix, iy; 
    ADC_SETUP_T  adc_hk;

    /* check to see if file exists */ 
    if (dmDatasetAccess(inp_p->adc_file, "R") == dmTRUE)
    { 
       short ii; 
       char  name[DS_SZ_PATHNAME];

       adc_hk.dataset = dmDatasetOpen(inp_p->adc_file);
       adc_hk.extension = dmBlockOpen(adc_hk.dataset, HPE_ADC_EXTNAME);
       adc_hk.num_cols = dmTableGetNoCols(adc_hk.extension);
       adc_hk.num_rows = dmTableGetNoRows(adc_hk.extension);
       adc_hk.row_check = dmSUCCESS;
  
       /* allocate buffers to store mapping/column data */
       if  (((adc_hk.mapping = 
            (short*) calloc(adc_hk.num_cols, sizeof(short))) != NULL) &&
            ((adc_hk.desc = (dmDescriptor**) calloc(adc_hk.num_cols,
             sizeof(dmDescriptor*))) != NULL))
       {
          for (ii = 0; ii < adc_hk.num_cols; ii++)
          {
             adc_hk.desc[ii] = 
                dmTableOpenColumnNo(adc_hk.extension, (ii+1));
             dmGetName(adc_hk.desc[ii], name, DS_SZ_KEYWORD);
             adc_hk.mapping[ii] = parse_adc_column(name);
          } /* end for */

          if (adc_hk.num_rows != (inp_p->x_taps + inp_p->y_taps))   
          {
             dsErrAdd(err_p, dsHPEADCROWCNTERR, Individual, Generic,
                      inp_p->adc_file); 
          } 
          
          /* add dependency check to ensure all data columns exist */ 
       
          /* load data from file */ 
          load_adc_entry(&adc_hk, inp_p, adc_x, adc_y, err_p);  

          /* check that all x taps read from adc correction table file */ 
          for (ix = 0; ix < inp_p->x_taps; ix++)
          {
             if (adc_x[ix].tap_num != ix)
             {
                dsErrAdd(err_p, dsHPEADCMISSROWERR, Individual, Generic,
                   inp_p->adc_file, 'u', ix);
                adc_x[ix].tap_num = ix; 
                adc_x[ix].p1 = adc_x[ix].p2 = adc_x[ix].p3 = 0; 
                adc_x[ix].q1 = adc_x[ix].q2 = adc_x[ix].q3 = 1; 
             }
          } 
          /* check that all y taps read from adc correction table file */ 
          for (iy = 0; iy < inp_p->y_taps; iy++)
          {
             if (adc_y[iy].tap_num != iy)
             {
                dsErrAdd(err_p, dsHPEADCMISSROWERR, Individual, Generic,
                   inp_p->adc_file, 'v', iy);
                adc_y[iy].tap_num = iy; 
                adc_y[iy].p1 = adc_y[iy].p2 = adc_y[iy].p3 = 0; 
                adc_y[iy].q1 = adc_y[iy].q2 = adc_y[iy].q3 = 1; 
             }
          } 
        
          /* close table */ 
          if (adc_hk.mapping != NULL)
          {
             free(adc_hk.mapping);
          }
          if (adc_hk.desc != NULL)
          {
             free(adc_hk.desc);
          }

          dmBlockClose(adc_hk.extension);
          dmDatasetClose(adc_hk.dataset);
       }
       else
       { 
          /* memory allocation failed  */
          dsErrAdd(err_p, dsHPEADCOPENMEMERR, Individual, Generic,
                   inp_p->adc_file); 
       } 
    }
    else
    {
       /* unable to access file */ 
       dsErrAdd(err_p, dsHPEADCOPENERR, Individual, Generic,
                inp_p->adc_file); 
    } 
}



/*H***********************************************************************

* DESCRIPTION:
 
  The routine parse_adc_column() is a private function called by 
  adc_table_load() to map the column names of the input adc correction
  file to index values so that the data can be properly loaded from
  the file regardless of the input file column orderings. The routine
  converts column names to uppercase before mapping them to allow
  case insensitivity on input file column names. If an unrecognized
  column name is encountered, it is mapped to the value ADC_UNKNOWN_COL
  (defined in adc_corr_defs.h). Unrecognized column names will be ignored
  when the data is read in.  

*H***********************************************************************/
short parse_adc_column(char* name) /* I - adc correction column name */
{
   short pos;
   short adc_mapping; 
   char name_uc[DS_SZ_PATHNAME];
 
   /* convert name to lowercase before comparing */
   if ((pos = strlen(name)) > DS_SZ_PATHNAME)
   {
      pos = DS_SZ_PATHNAME - 1;
   }
   name_uc[pos] = '\0';
   while (pos-- > 0)
   {
      name_uc[pos] = (char) toupper(name[pos]);
   }

   if ((strcmp(name_uc, ADC_AXIS_NAM)) == 0)
   {
      adc_mapping = ADC_AXIS_COL;
   }
   else if ((strcmp(name_uc, ADC_TAP_NAM)) == 0)
   {
      adc_mapping = ADC_TAP_COL;  
   }
   else if ((strcmp(name_uc, ADC_P1_NAM)) == 0)
   {
      adc_mapping = ADC_P1_COL;
   }
   else if ((strcmp(name_uc, ADC_Q1_NAM)) == 0)
   {
      adc_mapping = ADC_Q1_COL;
   }
   else if ((strcmp(name_uc, ADC_P2_NAM)) == 0)
   {
      adc_mapping = ADC_P2_COL;
   }
   else if ((strcmp(name_uc, ADC_Q2_NAM)) == 0)
   {
      adc_mapping = ADC_Q2_COL;
   }
   else if ((strcmp(name_uc, ADC_P3_NAM)) == 0)
   {
      adc_mapping = ADC_P3_COL;
   }
   else if ((strcmp(name_uc, ADC_Q3_NAM)) == 0)
   {
      adc_mapping = ADC_Q3_COL;
   }
   else
   {
      /* unknown column */
      adc_mapping = ADC_UNKNOWN_COL;
   }

   return (adc_mapping);
}
 


/*H***********************************************************************
 
* FUNCTION NAME: load_adc_entry()
 
* DESCRIPTION:
  This function takes in a mapping value (produced in parse_adc_column) 
  and reads the appropriate column from the current row of the adc
  correction file into the appropriate field of the adc correction data 
  structure. Unrecognized columns are ignored.
 
*H***********************************************************************/
 
void load_adc_entry (
   ADC_SETUP_P_T   hk_p,     /* I  - table pointer */
   INPUT_PARMS_P_T inp_p,    /* I  - input params  */
   ADC_CORR_P_T    adc_x,    /* O  - x tap values  */
   ADC_CORR_P_T    adc_y,    /* O  - y tap values  */
   dsErrList*      err_p)    /* O  - error list    */ 
{
   short mapping;
 
   hk_p->curr_row = 0; 

   while ((hk_p->curr_row < hk_p->num_rows) &&
       (hk_p->row_check != dmNOMOREROWS))
   {
      float p1=0.0, q1=0.0, p2=0.0, q2=0.0, p3=0.0, q3=0.0;
      short tap_num=0; 
      char axis[2];
       
      for (mapping = 0; mapping < hk_p->num_cols; mapping++)
      {
         switch (hk_p->mapping[mapping])
         {
            case ADC_AXIS_COL:
            {
                memset(axis, 0, 2);
                dmGetScalar_c(hk_p->desc[mapping], axis, 1); 
            }
            break;
 
            case ADC_TAP_COL:
               tap_num = dmGetScalar_s(hk_p->desc[mapping]);
            break;

            case ADC_P1_COL:
               p1 = dmGetScalar_f(hk_p->desc[mapping]); 
            break; 

            case ADC_Q1_COL:
               q1 = dmGetScalar_f(hk_p->desc[mapping]); 
            break; 

            case ADC_P2_COL:
               p2 = dmGetScalar_f(hk_p->desc[mapping]); 
            break; 

            case ADC_Q2_COL:
               q2 = dmGetScalar_f(hk_p->desc[mapping]); 
            break; 

            case ADC_P3_COL:
               p3 = dmGetScalar_f(hk_p->desc[mapping]);
            break;
 
            case ADC_Q3_COL:
               q3 = dmGetScalar_f(hk_p->desc[mapping]);
            break;

            default:
               /* do nothing */
            break;
         } /* end switch */
      } /* end for */
 
      switch (axis[0])
      {
         case 'u':   /* fallthrough intended */
         case 'U':   /* fallthrough intended */
            if ((tap_num > -1) && (tap_num < inp_p->x_taps))
            {
               adc_x[tap_num].tap_num = tap_num;
               adc_x[tap_num].p1 = p1;
               adc_x[tap_num].q1 = q1;
               adc_x[tap_num].p2 = p2;
               adc_x[tap_num].q2 = q2;
               adc_x[tap_num].p3 = p3;
               adc_x[tap_num].q3 = q3;
            }
         break;
         case 'v':   /* fallthrough intended */
         case 'V':   /* fallthrough intended */
            if ((tap_num > -1) && (tap_num < inp_p->y_taps))
            {
               adc_y[tap_num].tap_num = tap_num;
               adc_y[tap_num].p1 = p1;
               adc_y[tap_num].q1 = q1;
               adc_y[tap_num].p2 = p2;
               adc_y[tap_num].q2 = q2;
               adc_y[tap_num].p3 = p3;
               adc_y[tap_num].q3 = q3;
            }
         break;
         default:
           /* unknown/unexpected value- error */
           dsErrAdd(err_p, dsHPEADCLOADERR, Accumulation, Custom,
                    "ERROR: The axes values in %s must be either U or V.",
                    inp_p->adc_file); 
         break;
      }

      hk_p->row_check = dmTableNextRow(hk_p->extension);
      hk_p->curr_row++; 
   }
}


/*H***********************************************************************

* DESCRIPTION: 
 
  This module is called by hrc_process_events() to allocate the dynamic 
  memory needed to store the x and y axis adc correction tables. The size 
  of these tables varies with the type of system being used. 

  ie. in imaging mode the x adc correction table is 64 elements and 
                      the y adc correction table is 64 elements  

      in spectral mode the x adc correction table is 16 elements (TBR) and
                      the y adc correction table is 192 elements (TBR)    

      in spectral in imaging mode, the x adc correction table is 16 
                      elements (TBR) and the y adc correction table is 64
                      elements (TBR) 
* NOTES:

  The hrc system types for spectral and spectral with imaging need to be 
  firmly defined. Currently these are defined as "HRC-S" and "HRC-SI".  
 
  The routine returns a value of FALSE if no error was detected and TRUE if
  an error occurred. If errors are detected, they are added to the error list. 

  The routine assumes that the values for x_taps and y_taps in the
  input parameters data structure have already been correctly set. If not
  (as determined by checking for non-zero values), it will set the x_taps
  and y_taps in the same manner as the degap correction code.  

*H**************************************************************************/

boolean allocate_adc_table (
   INPUT_PARMS_P_T  inp_p,        /* I  - input parameters data structure  */
   ADC_CORR_P_T     *adc_x,       /* O  - array of x coord adc corrections */
   ADC_CORR_P_T     *adc_y,       /* O  - array of y coord adc corrections */
   dsErrList*       err_p)        /* O  - error list pointer               */
{
   boolean err_occurred = FALSE;

   /* set max_tap if needed (should be set by degap routines) */
   if ((inp_p->x_taps == 0) || (inp_p->y_taps == 0))  
   {
      if (inp_p->hrc_system == HRC_IMG_SYS)    /* imaging mode (hrc-i) */ 
      {
         inp_p->x_taps = HRC_I_X_TAPS;
         inp_p->y_taps = HRC_I_Y_TAPS;
         inp_p->min_tap[HDET_PLANE_X] = HRC_I_MIN_CRSU; 
         inp_p->max_tap[HDET_PLANE_X] = HRC_I_MAX_CRSU; 
         inp_p->min_tap[HDET_PLANE_Y] = HRC_I_MIN_CRSV; 
         inp_p->max_tap[HDET_PLANE_Y] = HRC_I_MAX_CRSV; 
      }
      else if (inp_p->hrc_system == HRC_SPC_SYS)  /*spectral mode */
      {
         inp_p->x_taps = HRC_S_X_TAPS;
         inp_p->y_taps = HRC_S_Y_TAPS;
         inp_p->hrc_system = HRC_SPC_SYS; 
         inp_p->min_tap[HDET_PLANE_X] = HRC_S_MIN_CRSU; 
         inp_p->max_tap[HDET_PLANE_X] = HRC_S_MAX_CRSU; 
         inp_p->min_tap[HDET_PLANE_Y] = HRC_S_MIN_CRSV; 
         inp_p->max_tap[HDET_PLANE_Y] = HRC_S_MAX_CRSV; 
      } 
      else if (inp_p->hrc_system == HRC_SPC_IMG_SYS) /*spectral imaging mode*/
      {
         inp_p->x_taps = HRC_SI_X_TAPS;
         inp_p->y_taps = HRC_SI_Y_TAPS;
         inp_p->min_tap[HDET_PLANE_X] = HRC_SI_MIN_CRSU; 
         inp_p->max_tap[HDET_PLANE_X] = HRC_SI_MAX_CRSU; 
         inp_p->min_tap[HDET_PLANE_Y] = HRC_SI_MIN_CRSV; 
         inp_p->max_tap[HDET_PLANE_Y] = HRC_SI_MAX_CRSV; 
      }
      else if (inp_p->hrc_system == HSI_IMG_SYS)  /* hsi (xrcf) */ 
      {
         inp_p->x_taps = HSI_X_TAPS;
         inp_p->y_taps = HSI_Y_TAPS;
         inp_p->min_tap[HDET_PLANE_X] = HSI_MIN_CRSU; 
         inp_p->max_tap[HDET_PLANE_X] = HSI_MAX_CRSU; 
         inp_p->min_tap[HDET_PLANE_Y] = HSI_MIN_CRSV; 
         inp_p->max_tap[HDET_PLANE_Y] = HSI_MAX_CRSV; 
      } 
      else
      {
         /* unknown hrc system mode specified */ 
         dsErrAdd(err_p, dsHPEADCBADSYSERR, Individual, Generic,
                  inp_p->instrume);
         err_occurred = TRUE;
      } 
   }

   if (!err_occurred) 
   {
      *adc_x = (ADC_CORR_P_T) calloc (inp_p->x_taps, sizeof(ADC_CORR_T));
      *adc_y = (ADC_CORR_P_T) calloc (inp_p->y_taps, sizeof(ADC_CORR_T));
   }
 
   if ((adc_x == NULL) || (adc_y == NULL))
   {
      /* calloc failed */ 
      dsErrAdd(err_p, dsHPEADCALLOCERR, Individual, Generic); 
      err_occurred = TRUE;
   }

   return (err_occurred); 
}



/*H***********************************************************************
 
* DESCRIPTION:
 
  The module deallocate_adc_table() is called by hrc_process_events to 
  free up the dynamic memory allocated for the x and y axis adc correction
  correction factor tables. After freeing the memory, it sets the pointers
  to NULL.  
 
*H***********************************************************************/

void deallocate_adc_table (
    ADC_CORR_P_T *adc_x,      /* I  - ptr to x adc correction table */
    ADC_CORR_P_T *adc_y )     /* I  - ptr to y adc correction table */
{
    if (*adc_x != NULL)
    {
       free (*adc_x);
       *adc_x = NULL; 
    }
    if (*adc_y != NULL)
    {
       free (*adc_y);
       *adc_y = NULL; 
    }
}



/*H***********************************************************************
 
* DESCRIPTION:

  The routine apply_adc_correction() is called by hrc_process_events to 
  perform the adc corrections on the raw amps. It takes in pointers to 
  the already loaded adc correction tables as well as a pointer to the 
  input parameters data structure and the current event record stucture.

  The tap positions of the current event are checked to see if they are
  within the accepted range (as specified in the input_parms_t data
  structure). The range is based upon the type of data being processed.
  If the tap positions fall within the accepted range, the adc
  correction values p1-p3 and q1-q3 are applied to the appropraite amps
  using the equation:

            AMPcorrected = Pn + Qn * AMPraw

* NOTES 

  The routine does not perform any error checking outside of the range 
  check. No error status is returned to the calling routine.  

*************************************************************************/

void apply_adc_correction(
    ADC_CORR_P_T  adc_x,      /* I  - ptr to x adc correction table     */
    ADC_CORR_P_T  adc_y,      /* I  - ptr to y adc correction table     */
    INPUT_PARMS_P_T inp_p,    /* I  - ptr to input parameters           */ 
    EVENT_REC_P_T evt_p)      /* I/O- event record structure pointer    */
{
   short cpx = evt_p->cp[HDET_PLANE_X],  /* crsu tap position */ 
         cpy = evt_p->cp[HDET_PLANE_Y];  /* crsv tap position */  

   if ((cpx >= inp_p->min_tap[HDET_PLANE_X]) && 
       (cpx <= inp_p->max_tap[HDET_PLANE_X]))
   {
      /* crsu within accepte4d range- apply adc corrections */ 
      evt_p->amps_dd[HDET_PLANE_X][HDET_1ST_AMP] = (adc_x[cpx].p1 +
         adc_x[cpx].q1 * evt_p->amps_dd[HDET_PLANE_X][HDET_1ST_AMP]);   
      evt_p->amps_dd[HDET_PLANE_X][HDET_2ND_AMP] = (adc_x[cpx].p2 +
         adc_x[cpx].q2 * evt_p->amps_dd[HDET_PLANE_X][HDET_2ND_AMP]);   
      evt_p->amps_dd[HDET_PLANE_X][HDET_3RD_AMP] = (adc_x[cpx].p3 +
         adc_x[cpx].q3 * evt_p->amps_dd[HDET_PLANE_X][HDET_3RD_AMP]);   
   }
   
   if ((cpy >= inp_p->min_tap[HDET_PLANE_Y]) && 
       (cpy <= inp_p->max_tap[HDET_PLANE_Y]))
   {
      /* crsv within accepte4d range- apply adc corrections */ 
      evt_p->amps_dd[HDET_PLANE_Y][HDET_1ST_AMP] = (adc_y[cpy].p1 +
         adc_y[cpy].q1 * evt_p->amps_dd[HDET_PLANE_Y][HDET_1ST_AMP]);   
      evt_p->amps_dd[HDET_PLANE_Y][HDET_2ND_AMP] = (adc_y[cpy].p2 +
         adc_y[cpy].q2 * evt_p->amps_dd[HDET_PLANE_Y][HDET_2ND_AMP]);   
      evt_p->amps_dd[HDET_PLANE_Y][HDET_3RD_AMP] = (adc_y[cpy].p3 +
         adc_y[cpy].q3 * evt_p->amps_dd[HDET_PLANE_Y][HDET_3RD_AMP]);   
   } 
} 


