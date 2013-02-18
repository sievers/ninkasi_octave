function[value]=calibrate_data_opts(tod,opts)
if (tod_has_calib_facs_c(tod))
  mdisp('using pre-stored calibration factors.');
  apply_calib_facs_c(tod);
  return
end

calib_facs=get_tod_calib_factors_opts(tod,opts);
mdisp('applying calibration');
apply_calib_facs_c(tod,calib_facs);
