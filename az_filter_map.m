function[map2]=az_filter_map(map,tods,invmap,myopts)
do_window=get_struct_mem(myopts,'window_symmetric',false);

map2=make_map_copy(map);
clear_map(map2);
for j=1:length(tods),
  mytod=tods(j);

  precalc_actpol_pointing_exact(mytod);
  convert_saved_pointing_to_pixellization(mytod,map);
  free_tod_pointing_saved(mytod);
  allocate_tod_storage(mytod);
  assign_tod_value(mytod,0.0);

  if (is_map_polarized(map))
    polmap2tod(map,mytod);
  else
    map2tod(map,mytod);
  end

  fit_hwp_az_poly_to_data(mytod,myopts);
  if (do_window)
    window_data(mytod);
    window_data(mytod);
  end

  
  if (is_map_polarized(map2))
    tod2polmap(mytod,map2);
  else
    tod2map(mytod,map2);
  end
  
  free_saved_pixellization(mytod);
  free_tod_storage(mytod);
end

mpi_reduce_map(map2);

if ~isempty(invmap)
  if (is_map_polarized(map2))
    apply_pol_precon_c(map2,invmap);
  else
    mm=skymap2octave(map2);
    wt=skymap2octave(invmap);
    mm=mm.*wt;
    octave2skymap(mm,map2);
  end
end

