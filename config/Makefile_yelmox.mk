#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

## EXTERNAL LIBRARIES #######################################

$(objdir)/basal_hydrology.o: $(libdir)/basal_hydrology.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/basal_hydro_simple.o: $(libdir)/basal_hydro_simple.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/geothermal.o: $(libdir)/geothermal.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/hyster.o: $(libdir)/hyster.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/latinhypercube.o: $(libdir)/latinhypercube.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/marine_shelf.o: $(libdir)/marine_shelf.f90 $(objdir)/pico.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/ismip6.o: $(libdir)/ismip6.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# $(objdir)/isostasy.o: $(libdir)/isostasy.f90 $(objdir)/nml.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

#$(objdir)/ncio.o: $(libdir)/ncio.f90
#	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

# $(objdir)/nml.o: $(libdir)/nml.f90
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/gaussian_filter.o: $(libdir)/gaussian_filter.f90
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/sediments.o: $(libdir)/sediments.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snapclim.o: $(libdir)/snapclim.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/stommel.o: $(libdir)/stommel.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/timeout.o: $(libdir)/timeout.f90 $(objdir)/nml.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/timer.o: $(libdir)/timer.f90
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/timestepping.o: $(libdir)/timestepping.f90 $(objdir)/nml.o $(objdir)/ncio.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/varslice.o: $(libdir)/varslice.f90 $(objdir)/nml.o $(objdir)/ncio.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# $(objdir)/xarray.o: $(libdir)/xarray.f90 $(objdir)/nml.o $(objdir)/ncio.o
# 	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# insol library 
$(objdir)/interp1D.o: $(libdir)/insol/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(libdir)/insol/insolation.f90 $(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<


# smbpal library
$(objdir)/smbpal_precision.o: $(libdir)/smbpal/smbpal_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_itm.o: $(libdir)/smbpal/smb_itm.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_pdd.o: $(libdir)/smbpal/smb_pdd.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(libdir)/smbpal/interp_time.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smbpal.o: $(libdir)/smbpal/smbpal.f90 $(objdir)/smbpal_precision.o $(objdir)/insolation.o  \
					$(objdir)/interp1D.o  $(objdir)/interp_time.o \
					$(objdir)/smb_pdd.o $(objdir)/smb_itm.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# pico library
$(objdir)/pico_geometry.o: $(libdir)/pico/pico_geometry.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/pico_physics.o: $(libdir)/pico/pico_physics.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/pico.o: $(libdir)/pico/pico.f90 $(objdir)/pico_geometry.o $(objdir)/pico_physics.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# General yelmox helper modules for different applications 
$(objdir)/yelmox_hysteresis_help.o: yelmox_hysteresis_help.f90 $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_YELMO) -c -o $@ $^

#############################################################
##							
## List of library files
##
#############################################################

yelmox_libs = 			$(objdir)/basal_hydrology.o \
						$(objdir)/basal_hydro_simple.o \
					    $(objdir)/geothermal.o \
					    $(objdir)/hyster.o \
					    $(objdir)/interp1D.o \
					    $(objdir)/insolation.o \
					    $(objdir)/ismip6.o \
					    $(objdir)/marine_shelf.o \
                        $(objdir)/pico_geometry.o \
                        $(objdir)/pico_physics.o \
                        $(objdir)/pico.o \
					    $(objdir)/sediments.o \
			 		    $(objdir)/smbpal_precision.o \
						$(objdir)/interp_time.o \
					    $(objdir)/smb_itm.o \
					    $(objdir)/smb_pdd.o \
					    $(objdir)/smbpal.o \
					    $(objdir)/snapclim.o

yelmox_help = 			$(objdir)/yelmox_hysteresis_help.o
