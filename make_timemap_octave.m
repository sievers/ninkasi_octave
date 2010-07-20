function[value]=make_timemap_octave(tods,map,mytimes)
%place the time difference relative to a standard in the data and project.

clear_map(map);

for j=1:length(tods), 
  ndet=get_tod_ndet(tods(j));
  n=get_tod_ndata(tods(j));
  t_start=get_tod_ctime_c(tods(j));
  dt=get_tod_dt(tods(j));

  tvec=dt*(0:(n-1) )'+t_start;
  tvec=tvec-mytimes(j);


  mydat=repmat(tvec,[1 ndet]);
  

  mytod=tods(j);
  allocate_tod_storage(mytod);
  push_tod_data(mydat,mytod);
  %assign_tod_value(mytod,1.0);
  %disp('assigned value');
  tod2map(mytod,map);
  %disp('projected');

  free_tod_storage(mytod);    
  %disp('freed');
end



  
