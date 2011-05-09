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


if isfield(a,'skymap')
  val=val+sum(sum(a.skymap.map.*b.skymap.map));
end

if isfield(a,'srccat')
  if iscell(a.srccat),
    for ss=1:numel(a.srccat),
      val=val+sum(a.srccat{ss}.amps.*b.srccat{ss}.amps);
    end
  else
    val=val+sum(a.srccat.amps.*b.srccat.amps);
  end
end

