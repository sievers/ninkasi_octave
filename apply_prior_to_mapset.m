function[mapset]=apply_prior_to_mapset(mapset,mapset_org)

if isfield(mapset,'skymap')
  if isfield(mapset.skymap,'has_prior')
    if (mapset.skymap.has_prior)
      disp('applying sky prior')
      mapset.skymap=apply_skymap_prior(mapset.skymap,mapset_org.skymap);
    end
  end
end

if isfield(mapset,'corrnoise')
  for j=1:length(mapset.corrnoise),
    if isfield(mapset.corrnoise(j),'prior')        
      %if (mapset.corrnoise(j).has_prior)
      %disp('applying correlated prior')
      mapset.corrnoise(j)=apply_corrnoise_prior(mapset.corrnoise(j));
      %end
    end
  end
end

function[map_in]=apply_skymap_prior(map_in,skymap)

map_in.map=map_in.map+sum(sum(skymap.map))*skymap.prior; 
%map_in.map=map_in.map+skymap.map*skymap.prior; 

return


function[corrnoise]=apply_corrnoise_prior(corrnoise)
mapft=fft(corrnoise.map);
for j=1:size(mapft,2),
  mapft(:,j)=mapft(:,j).*corrnoise.prior(:,j);
end

corrnoise.map=corrnoise.map+real(ifft(mapft));
return
