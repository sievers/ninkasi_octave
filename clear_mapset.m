function[mapset]=clear_mapset(mapset,copyptr)
%if copyptr comes in true, then we will allocate a new skymap in
%the pointer.  be careful doing this, as you can memory leak.

if ~exist('copyptr')
  copyptr=false;
end

if isfield(mapset,'skymap')
  mapset.skymap.map=0*mapset.skymap.map;
  if isfield(mapset.skymap,'mapptr')
    if (copyptr)
      mapset.skymap.mapptr=make_map_copy(mapset.skymap.mapptr);
    end
    clear_map(mapset.skymap.mapptr);
  end
end


if isfield(mapset,'cutvecs')
  if iscell(mapset.cutvecs)
    for j=1:numel(mapset.cutvecs),
      mapset.cutvecs(j)={0*mapset.cutvecs{j}};
    end
  else
    mapset.cutvecs=0*mapset.cutvecs;
  end
end


if isfield(mapset,'corrnoise'),
    for j=1:length(mapset.corrnoise),
        mapset.corrnoise(j).map=0*mapset.corrnoise(j).map;
    end
end

if isfield(mapset,'timestreams')
  for j=1:length(mapset.timestreams),
    mapset.timestreams(j).map=0*mapset.timestreams(j).map;
  end
end

if isfield(mapset,'srccat')
  if iscell(mapset.srccat),
    for ss=1:numel(mapset.srccat),
      mapset.srccat{ss}.amps=0*mapset.srccat{ss}.amps;
    end
  else
    mapset.srccat.amps=0*mapset.srccat.amps;
  end
end
