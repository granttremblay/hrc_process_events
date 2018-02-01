/*                                                                
**  Copyright (C) 1997-2007  Smithsonian Astrophysical Observatory 
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
 
* FILE NAME: bad_pixel_functions.c
 
* DEVELOPEMENT: tools
 
* DESCRIPTION: This file contains several functions used by 
  hrc_process_events to identify and flag hot spots (bad pixels).

  The modules in this file are:

     update_bad_pixel_list()
     cleanup_bad_pixel_data()
     check_for_bad_pixels()
     load_bad_pixel_files()
     open_bad_pixel_file()
     map_bad_pixel_column()
     load_bad_pixel_data()

  For details on the functionality of any of the above listed modules, 
  please see the 'description' comment preceding the specific module's 
  source code.
 
*H***********************************************************************/


#ifndef HRC_PROCESS_EVENTS_H
#include "hrc_process_events.h"
#define HRC_PROCESS_EVENTS_H
#endif 


/*************************************************************************

* FUNCTION NAME: update_bad_pixel_list()  

* DESCRIPTION:  The routine update_bad_pixel_list() is called on by
  load_bad_pixel_files() to insert hot spot (bad pixel) entries into the
  linked list pointed to by header pointer passed in as the second 
  function argument.

  The routine inserts elements into a linked list using the x coordinate 
  of the entry as a primary key. If more than one entry contains the same 
  x coordinate, the y coordinate is used as a second key. Duplicates 
  (items with the same x and y coordinates) are not maintained in the 
  list.   

* NOTES: 

  The routine frees the memory pointed to by the current hot pixel 
  pointer (first argument to function) if the linked list already 
  contains an entry with the same x and y coordinates. 

*************************************************************************/

void update_bad_pixel_list (
   BAD_PIX_P_T curr_p,  /* I  - current hot pixel    */
   BAD_PIX_P_T* head_p) /* I/O - hot pixel list head */
{
   int found = FALSE; 
   BAD_PIX_P_T prev_p = NULL;
   BAD_PIX_P_T temp_p = *head_p;

   /* iterate through the linked list to find insertion position */
   while (!found)
   {
      if (temp_p == NULL)
      { 
         /* at the end of the list */ 
         found = TRUE;
      }
      else if ((temp_p->x[0] < curr_p->x[0]) || 
               ((temp_p->x[0] == curr_p->x[0]) && 
                (temp_p->y[0] < curr_p->y[0]))) 
      { 
         /* not upto correct position in list, move to next list element */
         prev_p = temp_p;
         temp_p = temp_p->next; 
      } 
      else
      {
         /* have found the correct position in the list */ 
         found = TRUE;
      } 
   } 

   if (prev_p == NULL)
   {
      /* insert at head of list */
      if (temp_p != NULL)
      {
         if ((temp_p->x[0] == curr_p->x[0]) && (temp_p->x[1] == curr_p->x[1]) &&
             (temp_p->y[0] == curr_p->y[0]) && (temp_p->y[1] == curr_p->y[1])) 
         {
            /* entry is already in list */
            free (curr_p);
         }
         else
         {
            curr_p->next = *head_p; 
            *head_p = curr_p; 
         } 
      } 
      else
      {
         curr_p->next = NULL;
         *head_p = curr_p;
      } 
   }
   else
   {
      /* insert after temp element */
      if (temp_p != NULL)
      {
         if ((temp_p->x[0] == curr_p->x[0]) && (temp_p->x[1] == curr_p->x[1]) &&
             (temp_p->y[0] == curr_p->y[0]) && (temp_p->y[1] == curr_p->y[1]))
         {
            /* entry is already in list */
            free (curr_p);
         }
         else
         {
            curr_p->next = temp_p;
            prev_p->next = curr_p; 
         }
      } 
      else
      {
         curr_p->next = NULL;
         prev_p->next = curr_p; 
      } 
   } 
}




/*************************************************************************

* FUNCTION NAME: cleanup_bad_pixel_data() 
 
* DESCRIPTION:  The routine cleanup_bad_pixel_data() is called by 

  hrc_process_events to deallocate the memory that was used to store the 
  linked list of hot spots. The routine takes a pointer to the head of a 
  linked list and does not return a value. The routine iterates through 
  the linked list and frees all of the memory allocated to the elements 
  of the linked list.   

*************************************************************************/

void cleanup_bad_pixel_data(BAD_PIX_A_T bad_pix_list) /* I/O */  
{
   int rr;
   BAD_PIX_P_T current_p; 
 
   /* iterate through list to deallocate memory */
   for (rr = 4; rr--; )
   {
      while (bad_pix_list[rr] != NULL)
      {
         current_p = bad_pix_list[rr]; 
         bad_pix_list[rr] = bad_pix_list[rr]->next; 
         free(current_p);
      }
   }
}




/*************************************************************************
 
* FUNCTION NAMES: load_bad_pixel_files() 

* DESCRIPTION:  The routine load_bad_pixel_files() is called by 
  hrc_process_events to read in a stack of hot spot (bad pixel) files and 
  create a linked list of the coordinates. 
 
  The routine accepts a character string and a head pointer for the linked 
  list. A character string value of "NONE" or "none" indicates that no bad
  pixel files are to be loaded. Otherwise the routine calls on various hot 
  pixel routines to open the files, read in the data, and put the data 
  into the linked list. 

* NOTES: 

  The routine returns a value of TRUE if at least one bad pixel file was 
  successfully opened.  A value of FALSE indicates that no bad pixel
  files were opened or that the bad pixel stack was set to "NONE"/"none".  
 
*************************************************************************/
 
boolean load_bad_pixel_files(
   char*         badpix_file,   /* I - bad pixel file name   */
   BAD_PIX_P_T*  hotpix_list_p) /* I/O - hot pix list        */
{
   BAD_PIX_SETUP_T  bad_hk;  /* bad pixel house keeping data */
   char*       bad_p;
   boolean     opened = FALSE; /* TRUE if file successfully opened */
 
   bad_hk.stack = NULL;
   if ((ds_strcmp_cis(badpix_file, "none") == 0) ||
       (strcmp(badpix_file, "\0") == 0))
   {
      /* user has not supplied a bad pixel map */ 
      bad_p = NULL;
   }
   else
   {
      bad_p = badpix_file;
   }
   bad_hk.stack = stk_build(bad_p);
   
   /* use all bad pixel files in stack */ 
   while ((bad_hk.file = stk_read_next(bad_hk.stack)) != NULL)
   {
      bad_hk.curr_row = 0;        /* current row  read from badp tbl  */
 
      if (opened |= (open_bad_pixel_file(&bad_hk)))
      {
         load_bad_pixel_data(&bad_hk, hotpix_list_p);
 
         /* close bad pixel file */
         close_bad_pixel_file(&bad_hk); 
      }

      if (bad_hk.file != NULL)
      {
         free(bad_hk.file);
         bad_hk.file = NULL; 
      }
   }
 
   if (bad_hk.stack != NULL)
   {
      stk_close(bad_hk.stack);  
   }

   return (opened); 
}




/*************************************************************************
 
* FUNCTION NAME: open_bad_pixel_file()

* DESCRIPTION
 
  The routine open_bad_pixel_file is called from within a level one hrc 
  tool. Its purpose is to open up the specified bad pixel file table and 
  setup the necessary internal information (ie. column mappings via a call 
  to map_bad_pixel_columns().
 
  The routine returns a value of TRUE if the bad pixel file was opened 
  successfully and returns a value of FALSE if the file could not be 
  opened (or does not exist).
 
*************************************************************************/
 
boolean open_bad_pixel_file(
   BAD_PIX_SETUP_P_T  hk_p) /* I/O - bad pixel housekeeping info */ 
{
   boolean successful = TRUE;
 
   /* make sure a file name is provided before processing */
   if(((strcmp(hk_p->file, "NONE") != 0) && 
      (strcmp(hk_p->file, "none") != 0)) &&
      (strcmp(hk_p->file, "\0") != 0))
   {
      if (dmDatasetAccess(hk_p->file, "R") == dmTRUE)
      {
         int ii;

         /* open file */
         if ((hk_p->dataset = dmDatasetOpen(hk_p->file)) != NULL)
         {
            if ((hk_p->extension = dmBlockOpen(hk_p->dataset,"BADPIX")) != NULL)
            {
               hk_p->num_cols = dmTableGetNoCols(hk_p->extension);
               hk_p->num_rows = dmTableGetNoRows(hk_p->extension);
               hk_p->row_check = dmSUCCESS;
            } 
            else
            {
               successful = FALSE; 
            } 
         } 
         else
         {
            successful = FALSE; 
         } 
 
         /* allocate buffers to store mapping/column data */
         if  (successful && ((hk_p->mapping =
              (short*) calloc(hk_p->num_cols, sizeof(short))) != NULL) &&
              ((hk_p->desc = (dmDescriptor**) calloc(hk_p->num_cols,
               sizeof(dmDescriptor*))) != NULL))
         {
            char  name[DS_SZ_PATHNAME];

            for (ii = 0; ii < hk_p->num_cols; ii++)
            {
               hk_p->desc[ii] =
                  dmTableOpenColumnNo(hk_p->extension, (ii+1));
               dmGetName(hk_p->desc[ii], name, DS_SZ_KEYWORD);
               hk_p->mapping[ii] = map_bad_pixel_column(name);
            } /* end for */

         }  /* end if alloc succeeded */ 
      }
      else
      { 
         successful = FALSE;
      } 
   } 
 
   return (successful);
}



/*H********************************************************************
 
* FUNCTION NAME: map_bad_pixel_column()
 
* DESCRIPTION:
    This function takes in a string column name and returns a short 
    value which is used to map the data from the column into an 
    badpixel structure.
 
* NOTES:
    The mappings are done to allow for the easy introduction of 
    additional columns as the relevant input file columns may change.  
 
*H********************************************************************/
 
short map_bad_pixel_column (
   char name[])     /* I - number of columns in event */
{
    short mapping;
 
    if (ds_strcmp_cis(name, BPIX_CHIPX_COL) == 0)
    {
       mapping = BPIX_CHIPX;
    }
    else if (ds_strcmp_cis(name, BPIX_CHIPY_COL) == 0)
    {
       mapping = BPIX_CHIPY;
    }
    else if (ds_strcmp_cis(name, BPIX_CHIPID_COL) == 0)
    {
       mapping = BPIX_CHIPID;
    }
    else if (ds_strcmp_cis(name, BPIX_STATUS_COL) == 0)
    {
       mapping = BPIX_STATUS;
    }
    else
    {
       mapping = BPIX_UNKNOWN;
    }
 
    return (mapping);
}
 


/*H***********************************************************************
 
* FUNCTION NAME: load_bad_pixel_data()
 
* DESCRIPTION:
    This function takes in a mapping value (produced in the routine 
    map_bad_pixel_column) and reads the appropriate column from the 
    badpixel file to copy the data into the badpixel data structure.

* NOTES: 
    The mappings are done to allow for the easy introduction of additional
    columns as the relevant input file columns may change.  
 
*H***********************************************************************/
 
void load_bad_pixel_data (
    BAD_PIX_SETUP_P_T hk_p,   /* I - bad pixel file house keeping data */
    BAD_PIX_P_T*      hotpix_list_p) /* O - hot pixel list pointer     */
{
    BAD_PIX_P_T new_entry;
    short chip_id = 0; 
    short col;
 
    while ((hk_p->curr_row < hk_p->num_rows) &&
           (hk_p->row_check != dmNOMOREROWS))
    {
       new_entry = (BAD_PIX_P_T) malloc(sizeof(BAD_PIX_T));

       for (col = 0; col < hk_p->num_cols; col++)
       {
          switch (hk_p->mapping[col])
          {
             case BPIX_CHIPX:
                dmGetArray_l(hk_p->desc[col], new_entry->x, 2); 
             break;
       
             case BPIX_CHIPY:
                dmGetArray_l(hk_p->desc[col], new_entry->y, 2); 
             break;

             case BPIX_CHIPID:
                chip_id = dmGetScalar_s(hk_p->desc[col]);
             break; 

             case BPIX_STATUS:
             {
                unsigned char bit_val[1]  = {0};
                dmGetArray_bit(hk_p->desc[col], bit_val, 1);
                new_entry->status = bit_val[0];
             } 
             break; 
       
             default:
                /* do nothing */
             break;
          }
       }

       if ((chip_id > -1) && (chip_id < 4)) /* valid chip id? */ 
       {
          /* insert retrieved information into linked list */
          update_bad_pixel_list(new_entry, &hotpix_list_p[chip_id]);
       } 

       hk_p->row_check = dmTableNextRow(hk_p->extension);
       hk_p->curr_row++;
    }
} 

  
 

 

/*************************************************************************

* FUNCTION NAME: check_for_bad_pixels() 
 
* DESCRIPTION: The routine check_for_bad_pixels() is called to determine 
  if an event is located on a hot spot (bad pixel). The routine receives 
  a pointer to the head of an ordered linked list of hot spots and a 
  pointer to an event data structure. The routine iterates through the 
  linked list and checks the hot spot coordinates against the event's chip 
  coordinates. If the coordinates match a hot spot status bit is set in 
  the event data structure. The routine does not pass a return value. The 
  entire linked list is not traversed for each event, instead the 
  processing stops when it finds a match or, its x coordinate is greater 
  than the events x coordinate (or equivalent to the events x coord and 
  the lists y coord is greater than the events y coord).   
 
* NOTES: 

  The routine converts the chip coordinates of an event from doubles to 
  shorts via calls to pixlib. This is the same conversion which is done 
  to convert the raw coords from doubles to shorts when the events are 
  written to the output event file. 

  The routine expects the linked list to be ordered using the x 
  coordinate as the primary key and the y coordinate as a secondary key 
  for multiple hot spot instances on the same x coordinate. The list must 
  be in ascending order and not contain duplicate entries. Using the 
  load_bad_pixel_files routine to create the linked list ensures that 
  these conditions are met.
 
*************************************************************************/

void check_for_bad_pixels(
   BAD_PIX_A_T   bad_pixel_p, /* I - pointer to a list of hot pixels    */
   EVENT_REC_P_T evt_p)       /* I/O - event record structure pointer   */
{
   BAD_PIX_P_T  curr_p = NULL; 
   boolean      done = FALSE;
   short        chip_x, chip_y; 

   if ((evt_p->chipid > -1) && (evt_p->chipid < 4)) /* range 0-3 */
   {
      curr_p = bad_pixel_p[evt_p->chipid]; 
   } 
   else
   {
      /* warning- chipid out of range- no bad pixel checking */
   } 

   /* convert raw coords from doubles to shorts */ 
   pix_double_to_short(evt_p->chippos[HDET_PLANE_X], FALSE, &chip_x);
   pix_double_to_short(evt_p->chippos[HDET_PLANE_Y], FALSE, &chip_y);

   while (!done)
   {
      if (curr_p == NULL)
      {
         done = TRUE;
      }
      else if (curr_p->x[0] <= chip_x)
      {
         if (curr_p->x[1] < chip_x)
         {
            /* move to next element in hot pixel list */ 
            curr_p = curr_p->next; 
         } 
         else if (curr_p->y[0] <= chip_y) 
         {
            if (curr_p->y[1] < chip_y)
            {
               /* move to next element in hot pixel list */ 
               curr_p = curr_p->next;
            }
            else 
            {
               /* raw event is a hot pixel */
               evt_p->status |= HDET_HOT_SPOT_STS;
               done = TRUE; 
            }
         } 
         else 
         {
            /* move to next element in hot pixel list */ 
            curr_p = curr_p->next;
         } 
      } 
      else 
      { 
         /* since list is ordered - no need to check rest of list */
         done = TRUE;
      }
   } /* end while */ 
}

/*H***********************************************************************
 
* FUNCTION NAME: close_bad_pixel_file()
 
* DESCRIPTION:
 
  The routine close_bad_pixel_file() is called upon to clean up memory 
  and close datamodel file descripters used by the bad pixel routines.
 
*************************************************************************/
 
void close_bad_pixel_file(
   BAD_PIX_SETUP_P_T hk_p)   /* I/O - bad pixel file housekeeping info  */
{
   /* read all entries from file - perform clean up */
   if (hk_p->mapping != NULL)
   {
      free(hk_p->mapping);
   }
   if (hk_p->desc != NULL)
   {
      free(hk_p->desc);
   }
   if (hk_p->file != NULL)
   {
      free(hk_p->file);
      hk_p->file = NULL;
   }
 
   dmBlockClose(hk_p->extension);
   dmDatasetClose(hk_p->dataset);
}

