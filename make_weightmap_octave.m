function[value]=make_weightmap_octave(tods,map)
clear_map(map);

for j=1:length(tods),    
  mytod=tods(j);
  allocate_tod_storage(mytod);
  assign_tod_value(mytod,1.0);
  %disp('assigned value');
  tod2map(mytod,map);
  %disp('projected');

  free_tod_storage(mytod);    
  %disp('freed');
end



  
