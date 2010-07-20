function[mapset]=apply_precon(mapset,precon)
if isfield(precon,'skymap'),  %apply a Jacoby preconditioner to the sky
    mapset.skymap.map(precon.skymap.map>0)=mapset.skymap.map(precon.skymap.map>0)./precon.skymap.map(precon.skymap.map>0);    
    %mapset.skymap.map=mapset.skymap.map.*precon.skymap.map;
end
return

if isfield(precon,'corrnoise')
  for j=1:length(precon.corrnoise),
    if (isfield(precon(j).corrnoise,'prior'))
      mapset.corrnoise(j).map= fit_corrnoise_wpriors(precon.corrnoise(j).vecs,mapset.corrnoise(j).map,precon.corrnoise(j).prior);
    else
      mapset.corrnoise(j).map= fit_corrnoise_wpriors(precon.corrnoise(j).vecs,mapset.corrnoise(j).map);
    end
  end
end

      