function[val]=skymap2octave(map)
if strcmp(class(map),'int64')
  val=skymap2octave_c(map);
  return
end
if isstruct(map)
  
  if isfield(map,'partition');
    val=mpi_reduce_map(map);
  else
    val=map;
    val.map=skymap2octave_c(map.mapptr);
  end
end

