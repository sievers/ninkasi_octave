function[a]=add_mapset(a,b,fac)
if ~exist('fac')
    fac=1;
end

a.skymap.map=a.skymap.map-b.skymap.map;

for j=1:length(a.corrnoise),
    a.corrnoise(j).map=a.corrnoise(j).map-b.corrnoise(j).map;
end
