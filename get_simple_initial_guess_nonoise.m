function[mapset]=get_simple_initial_guess_nonoise(tods,mapset,myopts)
mapset=clear_mapset(mapset,'true');

for j=1:length(tods),
  mytod=tods(j);

  reload_process_data_opts(mytod,myopts);
  data=get_tod_data(mytod);
  free_tod_storage(mytod);

  if isfield(mapset,'corrnoise')
    vecs=mapset.corrnoise(j).vecs;
    ata=vecs*vecs';
    rhs=vecs*data';
    ts=(inv(ata)*rhs)';
    mapset.corrnoise(j).map=ts;
    data=data-ts*vecs;
  end
  if isfield(mapset,'timestreams')
     vecs=pull_bad_timestreams(mytod);
     ata=vecs'*vecs;
     rhs=vecs'*data;
     coeffs=inv(ata)*rhs;
     mapset.timestreams(j).map=coeffs;
     data=data-vecs*coeffs;
  end
  allocate_tod_storage(mytod);
  push_tod_data(data,mytod);
  clear data;
  tod2map(mytod,mapset.skymap.mapptr);
  
  free_tod_storage(mytod);
  
end

mm=mpi_allreduce(skymap2octave(mapset.skymap.mapptr));


wt=make_map_copy(mapset.skymap.mapptr);
make_weightmap_octave(tods,wt);
weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
weight_inv=1./weight;
weight_inv(weight==0)=0;
mm=mm.*weight_inv;

destroy_map(wt);
mapset.skymap.map=mm;

