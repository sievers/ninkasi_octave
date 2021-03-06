function[a]=add_mapset(a,b,fac)
if ~exist('fac')
    fac=1;
end

if isfield(a,'skymap')
  a.skymap.map=a.skymap.map+fac*b.skymap.map;
end

if isfield(a,'cutvecs')
  if iscell(a.cutvecs)
    for j=1:numel(a.cutvecs),
      a.cutvecs(j)={a.cutvecs{j}+fac*b.cutvecs{j}};
    end
  else
    a.cutvecs=a.cutvecs+b.cutvecs*fac;
  end
end

if isfield(a,'jumps')
  if iscell(a.jumps)
    for j=1:numel(a.jumps),
      a.jumps(j)={a.jumps{j}+fac*b.jumps{j}};
    end
  else
    a.jumps=a.jumps+b.jumps*fac;
  end
end


if isfield(a,'ground')
  a.ground.ground=a.ground.ground+fac*b.ground.ground;
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

