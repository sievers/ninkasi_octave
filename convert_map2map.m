function[wt]=convert_map2map(tods,map_in,map_out)

clear_map(map_out);
wt=make_map_copy(map_out);
clear_map(wt);

for j=1:length(tods),
  allocate_tod_storage(tods(j));
  map2tod(map_in,tods(j));
  tod2map(tods(j),map_out);
  assign_tod_value(tods(j),1.0);
  tod2map(tods(j),wt);  
  free_tod_storage(tods(j));
end

mm=mpi_allreduce(skymap2octave(map_out));
ww=mpi_allreduce(skymap2octave(wt));

mm(ww>0)=mm(ww>0)./ww(ww>0);
octave2skymap(mm,map_out);
clear mm
clear ww;




if nargout==0
  mdisp('getting rid of weight map');
  destroy_map(wt);
end

