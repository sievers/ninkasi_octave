function[mapset]=tod2mapset(mapset,tod,which_tod)
mapset.skymap.map=mapset.skymap.map+tod2map(tod.data,tod.ind,size(mapset.skymap.map,2),size(mapset.skymap.map,1));
%mapset.skymap.map=mapset.skymap.map+newmap;

if isfield(mapset,'corrnoise'),
    mapset.corrnoise(which_tod)=tod2corrnoise(tod.data,mapset.corrnoise(which_tod));
end
