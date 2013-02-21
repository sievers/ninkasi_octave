function[value]=destroy_tod(tod)
%try to clean up as much as possible from a TOD

free_tod_storage(tod);
free_tod_pointing_saved(tod);
free_tod_timevec_c(tod);
free_tod_altaz_c(tod);



free_tod_pointer_c(tod);  %better be the last call!