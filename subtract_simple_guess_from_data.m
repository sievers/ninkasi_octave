function[data_org]=subtract_simple_guess_from_data(mytod,mapset,myind)
data_org=get_tod_data(mytod);
assign_tod_value(mytod,1);
window_data(mytod);
window_vec=get_tod_data(mytod);
window_vec=window_vec(:,1);
free_tod_storage(mytod);
data=data_org;

if isfield(mapset,'timestreams'),
  vecs=pull_bad_timestreams(mytod);
  for j=1:size(vecs,2),
    vecs(:,j)=vecs(:,j).*window_vec;
  end

  rhs=data'*vecs;
  lhs=vecs'*vecs;
  fitp=rhs*inv(lhs);
  data=data-vecs*fitp';  
end
if isfield(mapset,'corrnoise'),
  vecs=mapset.corrnoise(myind).vecs;
  map=data*vecs';
  data=data-map*vecs;
end
push_tod_data(data,mytod);