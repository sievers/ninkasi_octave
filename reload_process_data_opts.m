function[value]=reload_process_data_opts(tods,myopts)

myopts=set_default_mapping_opts(myopts);
debutter=get_struct_mem(myopts,'debutter');
deconvolve_tau=get_struct_mem(myopts,'deconvolve_tau');
dedark=get_struct_mem(myopts,'dedark');
do_array_detrend=get_struct_mem(myopts,'detrend_array');
do_detrend=get_struct_mem(myopts,'detrend');
scale_fac=get_struct_mem(myopts,'scale_factor');
do_calib=get_struct_mem(myopts,'do_calib');


for j=1:length(tods)
  mytod=tods(j);

  allocate_tod_storage(mytod);
  read_tod_data(mytod);
  if (scale_fac~=1)
    multiply_all_data(mytod,scale_fac);
  end

  
  if do_calib
    calibrate_data_opts(mytod,myopts);
  end
  
  
  if (dedark)
    array_detrend(mytod);
    gapfill_data_c(mytod);
    subtract_darks_from_tod_eig (mytod);
  end
  
  gapfill_data_c(mytod);

  
  if (do_detrend|do_array_detrend)
    if (do_array_detrend)
      mdisp('array detrending');
      array_detrend(mytod);
    else
      mdisp('detrending');
      detrend_data_c_c(mytod);
    end
  end

  if (debutter)
    mdisp('debutterworthing');
    debutterworth_c(mytod);
  end
  if (deconvolve_tau)
    mdisp('deconvolving time constants.');
    deconvolve_tod_time_constants_c(mytod);
  end
  window_data(mytod);

end





