function[data]=mapset2tod(mapset,tod,which_tod) 
%pass in a mapset that continas a skymap and correlated noise maps, and get
%the expected tod.  Also need to specify which tod you want.  Need the tod 
%because that contains the pointing info.

skymap=mapset.skymap;
data=skymap2tod(tod.ind,skymap.map);

if isfield(mapset,'corrnoise')
    corrnoise=mapset.corrnoise(which_tod);
    data=data+corrnoise2tod(corrnoise);
end

