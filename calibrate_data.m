function[value]=calibrate_data(tod,varargin)

calib_facs=get_tod_calib_factors(tod,varargin{:});
apply_calib_facs_c(tod,calib_facs);
