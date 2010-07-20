function[x]=mapset_axpy(x,y,a,b)
if exist('b'),
  x.skymap.map=x.skymap.map*b;
  if isfield(x,'corrnoise'),
    for j=1:length(x.corrnoise),
      x.corrnoise(j).map=x.corrnoise(j).map*b;
    end
  end
end


if ~exist('a'),
  a=1;
end

x.skymap.map=x.skymap.map+a*y.skymap.map;
if isfield(x,'corrnoise'),
  for j=1:length(x.corrnoise),
    x.corrnoise(j).map=x.corrnoise(j).map+a*y.corrnoise(j).map;
  end
end
