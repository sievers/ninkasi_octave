function[val]=mapsetdotmapset(a,b)
val=0;
if isfield(a,'corrnoise'),
    for j=1:length(a.corrnoise),
        val=val+sum(sum(a.corrnoise(j).map.*b.corrnoise(j).map));
    end
    val=mpi_allreduce(val);
end

if isfield(a,'timestreams'),
    for j=1:length(a.timestreams),
        val=val+sum(sum(a.timestreams(j).map.*b.timestreams(j).map));
    end
    val=mpi_allreduce(val);
end


val=val+sum(sum(a.skymap.map.*b.skymap.map));
