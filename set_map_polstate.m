function[value]=set_map_polstate(map,pol)
max_npol=set_map_polstate_c;

%zero-pad in case we're short specified polarizations
if length(pol)<max_npol
  pol(max_npol)=0;
end
pol
set_map_polstate_c(map,pol);

