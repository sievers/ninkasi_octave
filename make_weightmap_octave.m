function[value]=make_weightmap_octave(tods,map,do_window)
clear_map(map);

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


  mdisp('assigned value');
  tod2map(mytod,map);
  mdisp('projected');

  free_tod_storage(mytod);    
  %disp('freed');
end



  
