function[value]=find_out_of_bounds_tods(tods,map)
if isstruct(map)
  if isfield(map,'mapptr')
    map=map.mapptr;
  else
    if isfield(map,'skymap')
      map=map.skymap.mapptr;
    end
  end
end
assert(class(map)=='int64');  %make sure we're in pointerland

ntod=numel(tods);
max_ntod=mpi_allreduce(ntod,'max');
mdisp(['biggest # of TODs is ' num2str(max_ntod)]);
myid=mpi_comm_rank+1;
for j=1:max_ntod,

  mdisp(['working on TOD ' num2str(j) ' of ' num2str(max_ntod)]);
  mpi_barrier;
  if j<=ntod,
    fid=fopen(['tod_check_node_' num2str(myid) '.txt'],'w');
    fprintf(fid,'%3d %s\n',j,get_tod_name(tods(j)));
    fclose(fid);
  end
  mpi_barrier;
  if j<=ntod,
    allocate_tod_storage(tods(j));
    tod2map(tods(j),map);  
    map2tod(map,tods(j));
    free_tod_storage(tods(j));
  end
  mpi_barrier;
  mdisp(['Finished.']);
end
