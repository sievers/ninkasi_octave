function[value]=calibrate_data_opts(tod,opts)

calib_facs=get_tod_calib_factors_opts(tod,opts);
mdisp('applying calibration');
apply_calib_facs_c(tod,calib_facs);
