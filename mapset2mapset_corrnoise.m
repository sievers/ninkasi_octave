function[new_mapset]=mapset2mapset_corrnoise(tods,mapset)

new_mapset=clear_mapset(mapset);
ntod=length(tods);
if isfield(mapset,'corrnoise'),
    assert(ntod==length(mapset.corrnoise));
end

%map_nx=size(skymap.map,1);
%map_ny=size(skymap.map,2);
assert(strcmp(class(mapset.skymap.map),'double'));

for j=1:ntod,
    tod=tods(j);
    assert(strcmp(class(tod.ind),'int32'));
    tod.data=mapset2tod(mapset,tod,j);
    if isfield(tod,'noise')
        tod=filter_data(tod);
    end
    new_mapset=tod2mapset(new_mapset,tod,j);   
end

