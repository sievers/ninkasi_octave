function[mapset]=apply_precon(mapset,precon)
if isfield(mapset,'skymap')
  if isfield(precon,'skymap'),  %apply a Jacoby preconditioner to the sky
    mapset.skymap.map(precon.skymap.map>0)=mapset.skymap.map(precon.skymap.map>0)./precon.skymap.map(precon.skymap.map>0);    
    %mapset.skymap.map=mapset.skymap.map.*precon.skymap.map;
  end
end


%since we often cut windows, then use a preconditioner for the cut section.
if isfield(mapset,'cutvecs')
  if isfield(precon,'cutvecs')
    %mdisp('applying cutvecs preconditioner')
    if iscell(mapset.cutvecs)
      for j=1:numel(mapset.cutvecs),
        mapset.cutvecs(j)={mapset.cutvecs{j}./precon.cutvecs{j}};
      end
    else
      mapset.cutvecs=mapset.cutvecs./precon.cutvecs;
    end
  end
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

      