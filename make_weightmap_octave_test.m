function[value]=make_weightmap_octave_test(tods,map,do_window,varargin)

do_new_pointing=get_keyval_default('do_new_pointing',false,varargin{:})
do_reduce=get_keyval_default('mpi_reduce',false,varargin{:});

clear_map(map);

if isempty('do_window')
  clear do_window;
end


if ~exist('do_window')
  do_window=true;
end

for j=1:length(tods),    
  mytod=tods(j);

  if (do_new_pointing)
    mdisp('doing pointing')
    precalc_actpol_pointing_exact(mytod);
    convert_saved_pointing_to_pixellization(mytod,map);
    free_tod_pointing_saved(mytod);
  end
  


  allocate_tod_storage(mytod);
  assign_tod_value(mytod,1.0);
  if do_window,
    window_data(mytod);
    window_data(mytod);
  end

  
  
  mdisp('assigned value');
  if (is_map_polarized(map))
    tod2polmap(mytod,map);
  else
    tod2map(mytod,map);
  end
  mdisp('projected');
  
  if (do_new_pointing)
    %free_tod_pointing_saved(mytod);
    free_saved_pixellization(mytod);
    
  end


  free_tod_storage(mytod);    
  %disp('freed');
end

if (do_reduce)
  mdisp('reducing')
  weight=skymap2octave(map);
  weight=mpi_allreduce(weight);
  octave2skymap(weight,map);
  clear weight;
end

  
