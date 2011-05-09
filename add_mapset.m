function[a]=add_mapset(a,b,fac)
if ~exist('fac')
    fac=1;
end

if isfield(a,'skymap')
  a.skymap.map=a.skymap.map+fac*b.skymap.map;
end


if (isfield(a,'corrnoise')),
    for j=1:length(a.corrnoise),
        a.corrnoise(j).map=a.corrnoise(j).map+fac*b.corrnoise(j).map;
    end
end


if (isfield(a,'timestreams')),
    for j=1:length(a.timestreams),
        a.timestreams(j).map=a.timestreams(j).map+fac*b.timestreams(j).map;
    end
end

if isfield(a,'srccat'),
  a.srccat.amps=a.srccat.amps+fac*b.srccat.amps;
end

