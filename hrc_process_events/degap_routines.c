/*                                                                
**  Copyright (C) 1996-2007  Smithsonian Astrophysical Observatory 
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
 
* FILE NAME: degap_routines.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: 

  The file degap_routines.c contains the following modules used by 
  hrc_process_events() in dealing with degap tables:

        allocate_degap_table()  
        deallocate_degap_table()  
        degap_table_load()  
        parse_degap_column() 
 
* NOTES:

  wmclaugh@cfa	June 11, 1996  First Version.

* REVISION HISTORY:
 
        Ref. No.        Date
        --------        ----
        1.1             11 Jun 1996
*H**************************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#endif 

/*H**************************************************************************

* DESCRIPTION: 

  The module degap_table_load() is used by hrc_process_events to populate the 2
  degapping tables used by the finepos algorithm for degap corrections. The two 
  tables (one for the x axis and one for the y axis) may be populated in any of 
  the following three methods as specified by the degap parameter...
 
     NONE,none,"" - All linear correction factors are set to one and all
                    quadratic correction factors are set to zero
 
     COEFF        - The four correction factors specified as parameters (cfx1,
                    cfx2, cfy1, and cfy2) are used to populate the tables
 
     <filename>   - The correction factors are read from a table file
 
 
* NOTES:
 
  The routine does not return any status but errors are added onto the 
  error list passed in. 
 
*H**************************************************************************/

void degap_table_load(
    INPUT_PARMS_P_T  inp_p,     /* I   pointer to input parameters  */
    DEGAP_P_T     degap_x,      /* O   array of x coord degaps */ 
    DEGAP_P_T     degap_y,      /* O   array of y coord degaps */
    dsErrList*    err_p)        /* O   error list              */
{
    double lax, lbx, rax, rbx, lay, lby, ray, rby;
    short ix, iy; 
    boolean use_coeff; 

    if ((use_coeff = (ds_strcmp_cis(inp_p->degap_file, "COEFF") == 0)) ||  
        ((ds_strcmp_cis(inp_p->degap_file, "NONE") == 0) ||
         (strcmp(inp_p->degap_file, "\0") == 0)))
    {
       if (use_coeff)   /* use correction factors from parameter table */
       {
          lax = inp_p->cf[HDET_PLANE_X][HDET_1ST_ORD_CF];
          lbx = inp_p->cf[HDET_PLANE_X][HDET_2ND_ORD_CF];
          rax = inp_p->cf[HDET_PLANE_X][HDET_1ST_ORD_CF];
          rbx = inp_p->cf[HDET_PLANE_X][HDET_2ND_ORD_CF];
          lay = inp_p->cf[HDET_PLANE_Y][HDET_1ST_ORD_CF];
          lby = inp_p->cf[HDET_PLANE_Y][HDET_2ND_ORD_CF];
          ray = inp_p->cf[HDET_PLANE_Y][HDET_1ST_ORD_CF];
          rby = inp_p->cf[HDET_PLANE_Y][HDET_2ND_ORD_CF];
       }
       else /* use hard wired "NO DEGAP" values */ 
       {
          rax = lax = lay = ray = 1.0;
          rbx = lbx = lby = rby = 0.0;
       }

       /* populate x degap table */ 
       for (ix = 0; ix < inp_p->x_taps; ix++)
       {
          degap_x[ix].tap_num = ix;
          degap_x[ix].la = lax;
          degap_x[ix].lb = lbx;
          degap_x[ix].ra = rax;
          degap_x[ix].rb = rbx;
       }
       /* populate y degap table */ 
       for (iy = 0; iy < inp_p->y_taps; iy++)
       {
          degap_y[iy].tap_num = iy;
          degap_y[iy].la = lay;
          degap_y[iy].lb = lby;
          degap_y[iy].ra = ray;
          degap_y[iy].rb = rby;
       }
    }
    else       /* load data from file */ 
    {
       DEGAP_SETUP_T  dgp_hk;

       /* check to see if file exists */ 
       if (dmDatasetAccess(inp_p->degap_file, "R") == dmTRUE)
       { 
          short ii; 
          char  name[DS_SZ_PATHNAME];

          dgp_hk.dataset = dmDatasetOpen(inp_p->degap_file);
          dgp_hk.extension = dmBlockOpen(dgp_hk.dataset, HPE_DEGAP_EXTNAME);
          dgp_hk.num_cols = dmTableGetNoCols(dgp_hk.extension);
          dgp_hk.num_rows = dmTableGetNoRows(dgp_hk.extension);
          dgp_hk.row_check = dmSUCCESS;
  
          /* allocate buffers to store mapping/column data */
          if  (((dgp_hk.mapping = 
               (short*) calloc(dgp_hk.num_cols, sizeof(short))) != NULL) &&
               ((dgp_hk.desc = (dmDescriptor**) calloc(dgp_hk.num_cols,
                sizeof(dmDescriptor*))) != NULL))
          {
             for (ii = 0; ii < dgp_hk.num_cols; ii++)
             {
                dgp_hk.desc[ii] = 
                   dmTableOpenColumnNo(dgp_hk.extension, (ii+1));
                dmGetName(dgp_hk.desc[ii], name, DS_SZ_KEYWORD);
                dgp_hk.mapping[ii] = parse_degap_column(name);
             } /* end for */
 
             if (dgp_hk.num_rows != (inp_p->x_taps + inp_p->y_taps))   
             {
                dsErrAdd(err_p, dsHPEADCROWCNTERR, Individual, Generic,
                         inp_p->degap_file);
             } 
             else 
             {
                /* dependency check to ensure all data columns exist */ 

                load_degap_entry(&dgp_hk, inp_p, degap_x, degap_y, err_p);  

                /* check that all x taps read from degap table file */ 
                for (ix = 0; ix < inp_p->x_taps; ix++)
                {
                   if (degap_x[ix].tap_num != ix)
                   {
                      dsErrAdd(err_p, dsHPEADCMISSROWERR, Individual, Custom,
                        "WARNING: %s does not contain values for %c tap %hd.",
                        inp_p->degap_file, 'u', ix);
                   }
                } 
                /* check that all y taps read from degap table file */ 
                for (iy = 0; iy < inp_p->y_taps; iy++)
                {
                   if (degap_y[iy].tap_num != iy)
                   {
                      dsErrAdd(err_p, dsHPEADCMISSROWERR, Individual, Custom,
                        "WARNING: %s does not contain values for %c tap %hd.",
                        inp_p->degap_file, 'v', iy);
                   }
                } 
             } 
             /* close table */ 
             if (dgp_hk.mapping != NULL)
             {
                free(dgp_hk.mapping);
             }
             if (dgp_hk.desc != NULL)
             {
                free(dgp_hk.desc);
             }
 
             dmBlockClose(dgp_hk.extension);
             dmDatasetClose(dgp_hk.dataset);
          }
          else
          { 
             dsErrAdd(err_p, dsHPEDEGAPLOADERR, Individual, Generic); 
          } 
       }
       else
       {
          /* unable to open */
          dsErrAdd(err_p, dsHPEADCOPENMEMERR, Individual, Generic,
                   inp_p->degap_file);
       } 
    }
}

short parse_degap_column(char* name) /* I - degap column name */
{
   short pos;
   short degap_mapping; 
   char name_lc[DS_SZ_PATHNAME];
 
   /* convert name to lowercase before comparing */
   if ((pos = strlen(name)) > DS_SZ_PATHNAME)
   {
      pos = DS_SZ_PATHNAME - 1;
   }
   name_lc[pos] = '\0';
   while (pos-- > 0)
   {
      name_lc[pos] = (char) tolower(name[pos]);
   }

  if ((strcmp(name_lc, DGP_AXIS_NAM)) == 0)
  {
     degap_mapping = DGP_AXIS_COL;
  }
  else if ((strcmp(name_lc, DGP_TAP_NAM)) == 0)
  {
     degap_mapping = DGP_TAP_COL;  
  }
  else if ((strcmp(name_lc, DGP_LA_NAM)) == 0)
  {
     degap_mapping = DGP_LA_COL;
  }
  else if ((strcmp(name_lc, DGP_LB_NAM)) == 0)
  {
     degap_mapping = DGP_LB_COL;
  }
  else if ((strcmp(name_lc, DGP_RA_NAM)) == 0)
  {
     degap_mapping = DGP_RA_COL;
  }
  else if ((strcmp(name_lc, DGP_RB_NAM)) == 0)
  {
     degap_mapping = DGP_RB_COL;
  }
  else
  {
     /* unknown column */
     degap_mapping = HDET_UNKNOWN_FIELD;
  }

  return (degap_mapping);
}
 


/*H**************************************************************************
 
* FUNCTION NAME: load_degap_entry()
 
* DESCRIPTION:
    This function takes in a mapping value (produced in parse_degap_column) and
    reads the appropriate column from the current row of the degap file into
    the appropriate field of the degap data structure.
 
*H**************************************************************************/
 
void load_degap_entry (DEGAP_SETUP_P_T hk_p, /* table pointer */
                       INPUT_PARMS_P_T inp_p, /* input params */
                       DEGAP_P_T degap_x, /* O- x tap values  */
                       DEGAP_P_T degap_y, /* O- y tap values  */
                       dsErrList* err_p)  /* O- error list    */
{
   short mapping;
 
   hk_p->curr_row = 0; 

   while ((hk_p->curr_row < hk_p->num_rows) &&
       (hk_p->row_check != dmNOMOREROWS))
   {
      double la, lb, ra, rb; 
     short tap_num; 
      char axis[2];
       
      for (mapping = 0; mapping < hk_p->num_cols; mapping++)
      {
         switch (hk_p->mapping[mapping])
         {
            case DGP_AXIS_COL:
                dmGetScalar_c(hk_p->desc[mapping], axis, 1); 
 
            break;
 
            case DGP_TAP_COL:
               tap_num = dmGetScalar_s(hk_p->desc[mapping]);
            break;

            case DGP_LA_COL:
               la = dmGetScalar_d(hk_p->desc[mapping]); 
            break; 

            case DGP_LB_COL:
               lb = dmGetScalar_d(hk_p->desc[mapping]); 
            break; 

            case DGP_RA_COL:
               ra = dmGetScalar_d(hk_p->desc[mapping]); 
            break; 

            case DGP_RB_COL:
               rb = dmGetScalar_d(hk_p->desc[mapping]); 
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
               degap_x[tap_num].tap_num = tap_num;
               degap_x[tap_num].la = la;
               degap_x[tap_num].lb = lb;
               degap_x[tap_num].ra = ra;
               degap_x[tap_num].rb = rb;
            }
         break;
         case 'v':   /* fallthrough intended */
         case 'V':   /* fallthrough intended */
            if ((tap_num > -1) && (tap_num < inp_p->y_taps))
            {
               degap_y[tap_num].tap_num = tap_num;
               degap_y[tap_num].la = la;
               degap_y[tap_num].lb = lb;
               degap_y[tap_num].ra = ra;
               degap_y[tap_num].rb = rb;
            }
         break;
         default:
           /* unknown/unexpected value- error */
           dsErrAdd(err_p, dsHPEDEGAPLOADERR, Accumulation, Custom,
                    "ERROR: The axes values in %s must be either U or V.",
                    inp_p->degap_file);
         break;
      }

      hk_p->row_check = dmTableNextRow(hk_p->extension);
      hk_p->curr_row++; 
   }
}




/*H**************************************************************************

* DESCRIPTION: 
 
  This module is called by hrc_process_events() to allocate the dynamic memory
  needed to store the x and y axis degap tables. The size of these tables
  varies with the type of system being used. 

     ie. in imaging mode  the x degap table is 64 elements and 
                          the y degap table is 64 elements  

         in spectral mode the x degap table is 16 elements (TBR)  and
                          the y degap table is 192 elements (TBR)    

         in spectral in imaging mode, the x degap table is 16 elements (TBR) 
                          and the y degap table is 64 elements (TBR) 
 
 
* NOTES:

  The hrc system types for spectral and spectral with imaging need to be 
  firmly defined. Currently these are defined as "HRC-S" and "HRC-SI".  
 
  The routine returns a value of TRUE if an error occurred during
  execution. If no errors were detected a value of FALSE is returned.

*H**************************************************************************/

boolean allocate_degap_table (
    INPUT_PARMS_P_T inp_p,       /* I  - input params */
    DEGAP_P_T     *degap_x,      /* O   array of x coord degaps */
    DEGAP_P_T     *degap_y)      /* O   array of y coord degaps */
{
    boolean err_occurred = FALSE;

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
    else if (inp_p->hrc_system == HRC_SPC_IMG_SYS)  /*spectral imaging mode */
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
       err_occurred = TRUE; 
    }

    if (!err_occurred) 
    {
       *degap_x = (DEGAP_P_T) calloc (inp_p->x_taps, sizeof(DEGAP_T));
       *degap_y = (DEGAP_P_T) calloc (inp_p->y_taps, sizeof(DEGAP_T));
    }
 
    if ((degap_x == NULL) || (degap_y == NULL))
    {
       /* calloc failed */ 
       err_occurred = TRUE; 
    }

    return (err_occurred);
}


/*H**************************************************************************
 
* DESCRIPTION:
 
  The module deallocate_degap_table() is called by hrc_process_events to free up
  the dynamic memory allocated for the x and y axis degap correction
  factor tables. 
 
*H**************************************************************************/

void deallocate_degap_table (
    DEGAP_P_T *degap_x,      /* I  - ptr to x degapping table */
    DEGAP_P_T *degap_y )     /* I  - ptr to y degapping table */
{
    if (*degap_x != NULL)
    {
       free (*degap_x);
    }
    if (*degap_y != NULL)
    {
       free (*degap_y);
    }
}


