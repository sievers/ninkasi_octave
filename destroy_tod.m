function[value]=destroy_tod(tod)
%try to clean up as much as possible from a TOD

destroy_actpol_pointing(tod);
destroy_actpol_pointing_fit(tod);
free_tod_storage(tod);
free_tod_pointing_saved(tod);
free_tod_timevec_c(tod);
free_tod_altaz_c(tod);
free_tod_hwp_c(tod);
free_tod_jumps(tod);
%some more things to clean up - JLS 23/1/14
free_saved_pixellization(tod);
free_tod_data_saved(tod);

free_tod_pointer_c(tod);  %better be the last call!