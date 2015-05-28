function[map]=make_weightmap_octave(tods,map,do_window,varargin)

do_new_pointing=get_keyval_default('do_new_pointing',false,varargin{:})
do_reduce=get_keyval_default('mpi_reduce',false,varargin{:});
free_2gamma=get_keyval_default('free_2gamma',true,varargin{:});



if isstruct(map)
  if isfield(map,'partition')
    clear_map(map.mapptr);
    make_weightmap_octave(tods,map.mapptr,do_window,'mpi_reduce',false,varargin{:});
    if do_reduce
      map=mpi_reduce_map(map);
    end
    return;
  end
end



clear_map(map);

if isempty('do_window')
  clear do_window;
end


if ~exist('do_window')
  do_window=true;
end

for j=1:length(tods),    
  mytod=tods(j);
  allocate_tod_storage(mytod);
  assign_tod_value(mytod,1.0);
  if do_window,
    window_data(mytod);
    window_data(mytod);
  end

  
  if (do_new_pointing)
    mdisp('doing pointing')
    precalc_actpol_pointing_exact(mytod);
    convert_saved_pointing_to_pixellization(mytod,map);
    free_tod_pointing_saved(mytod,free_2gamma);
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

  
