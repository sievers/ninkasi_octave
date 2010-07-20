function[mapset]=create_initial_mapset(tods,mapset)
mapset=clear_mapset(mapset);
if isfield(mapset,'corrnoise'),
    assert(length(tods)==length(mapset.corrnoise));
end

for j=1:length(tods),
    tod=tods(j);
    if isfield(tod,'noise'),
        tod=filter_data(tod);
    end
    mapset=tod2mapset(mapset,tod,j);
end
