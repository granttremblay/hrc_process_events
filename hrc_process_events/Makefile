#-----------------------------------------------------------------------
#
# 		Makefile for the hrc_process_events program C verion
#
#-----------------------------------------------------------------------

MK_TOP = ../../..
include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

LOCAL_LIBS        = -L$(SCI_ROOT)/lib -ll1_asp -L../l1_hrc -ll1_hrc 
LOCAL_INC         = -I../l1_hrc

EXEC              = hrc_process_events
PAR_FILES         = hrc_process_events.par
XML_FILES         = hrc_process_events.xml

SRCS	= amp_sf_cor_functions.c \
          tap_ring_functions.c \
          calculate_pi_hrc.c \
          hpe_gain.c \
	  hpe_setup_degap_file.c \
          adc_corr_routines.c \
          badpixel_functions.c \
          coordinate_transforms.c \
          process_warnings.c \
          dependency_check_hrc.c \
          hrc_process_events.c \
          hrc_process_configure_pixlib.c \
          hrc_process_events_functions.c \
          hrc_process_read_obsfile.c \
          hrc_process_set_instrume.c \
          hrc_process_setup_input_file.c \
          hrc_process_setup_output_file.c \
          hrc_process_setup_logfile.c \
          hrc_process_time_check.c \
          hrc_setup_columns.c \
          load_event_data.c \
          parse_hrc_evt_columns.c \
          ratio_checks_hrc.c \
          sum_phas_hrc.c \
          t_hrc_process_events.c \
          write_hrc_events.c \
          write_instrume_params.c \
	  adc_filter_routines.c \
	  hpe_setup_calibration.c


OBJS	= $(SRCS:.c=.o)

MAKETEST_SCRIPT   = hrc_process_events.t


include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------
$(EXEC):$(OBJS) 
	$(LINK)
	@echo

announce1:
	@echo "   /---------------------------------------------------------\ "
	@echo "   |            Building hrc_process_events program          | "
	@echo "   \---------------------------------------------------------/ "

