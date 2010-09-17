function[mapset]=apply_prior_to_mapset_mustang(mapset,mapset_org)

if isfield(mapset,'skymap')
  if isfield(mapset.skymap,'has_prior')
    if (mapset.skymap.has_prior)
      disp('applying sky prior')
      mapset=apply_masked_skymap_prior(mapset,mapset_org);
    end
  end
end

if isfield(mapset,'corrnoise')
  for j=1:length(mapset.corrnoise),
    if isfield(mapset.corrnoise(j),'prior')        
      mapset.corrnoise(j)=apply_corrnoise_prior(mapset.corrnoise(j));
    end
  end
end



function[map_in]=apply_skymap_prior(map_in,skymap)
map_in.map=map_in.map+sum(sum(skymap.map))*skymap.prior; 
%map_in.map=map_in.map+skymap.map*skymap.prior; 
return




function[corrnoise]=apply_corrnoise_prior(corrnoise)

[mapft,x]=get_corrnoise_ft_nu(corrnoise.map,'do_abs',false);

for j=1:size(mapft,2),
  x(1)=x(2);  %gyrations to avoid an undefined value in the zero-frequency mode
  vec=corrnoise.prior(1,j)+corrnoise.prior(3,j)*(x.^corrnoise.prior(2,j)); %we're assuming a power law fit here
  vec=1./vec.^2; %should be one, if not, it's for testing.
  vec(1)=0;
  mapft(:,j)=mapft(:,j).*vec;
end
n=size(corrnoise.map,1);
corrnoise.map=corrnoise.map+fft_c2r(mapft,iseven(n))*n;

return
