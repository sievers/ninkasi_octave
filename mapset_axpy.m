function[x]=mapset_axpy(x,y,a,b)

if exist('b'),
  x.skymap.map=x.skymap.map*b;
  if isfield(x,'corrnoise'),
    for j=1:length(x.corrnoise),
      x.corrnoise(j).map=x.corrnoise(j).map*b;
    end
  end
else
  b=1;
end


if ~exist('a'),
  a=1;
end

x.skymap.map=x.skymap.map+a*y.skymap.map;
if isfield(x,'ground')
  x.ground.ground=b*x.ground.ground+a*y.ground.ground;
end


if isfield(x,'corrnoise'),
  for j=1:length(x.corrnoise),
    x.corrnoise(j).map=x.corrnoise(j).map+a*y.corrnoise(j).map;
  end
end

if isfield(x,'srccat')
  if iscell(x.srccat),
    for ss=1:numel(x.srccat),
      x.srccat{ss}.amps=x.srccat{ss}.amps+a*y.srccat{ss}.amps;
    end
  else
    x.srccat.amps=x.srccat.amps+a*y.srccat.amps;
  end
end
