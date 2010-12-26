function[mapset]=apply_masked_skymap_prior(mapset,mapset_org)

if isfield(mapset,'skymap')
  if isfield(mapset.skymap,'has_prior')
    if (mapset.skymap.has_prior)
      mapset.skymap=apply_skymap_prior(mapset.skymap,mapset_org.skymap);
    end
  end
end

function[map_in]=apply_skymap_prior(map_in,skymap)

%disp(sum(sum(skymap.map.*skymap.prior_mask)))
map_in.map=map_in.map-sum(sum(skymap.map.*skymap.prior_mask));

return

