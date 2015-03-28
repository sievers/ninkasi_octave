function[val]=mapsetdotmapset(a,b)
val=0;
if isfield(a,'corrnoise'),
    for j=1:length(a.corrnoise),
        val=val+sum(sum(a.corrnoise(j).map.*b.corrnoise(j).map));
    end
    val=mpi_allreduce(val);
end

if isfield(a,'timestreams'),
  tot=0;
  for j=1:length(a.timestreams),
    tot=tot+sum(sum(a.timestreams(j).map.*b.timestreams(j).map));
  end
  val=val+mpi_allreduce(tot);
end


if isfield(a,'skymap')
  tmp=sum(sum(sum(a.skymap.map.*b.skymap.map)));
  if isfield(a.skymap.map,'partition')
    tmp=mpi_allreduce(tmp);
  end
  val=val+tmp;
end

if isfield(a,'ground')
  %val=val+map_times_map(a.ground.groundptr,b.ground.groundptr);
  myval=sum(sum(sum(a.ground.ground.*b.ground.ground)));

  val=val+mpi_allreduce(myval);

end



if isfield(a,'cutvecs')
  tot=0;
  if iscell(a.cutvecs)
    for j=1:numel(a.cutvecs),
      tot=tot+sum( (a.cutvecs{j}).*(b.cutvecs{j}));
    end
  else
    tot=tot+sum(a.cutvecs.*b.cutvecs);
  end
  val=val+mpi_allreduce(tot);
end


if isfield(a,'srccat')
  tot=0;
  if iscell(a.srccat),
    for ss=1:numel(a.srccat),
      tot=tot+sum(a.srccat{ss}.amps.*b.srccat{ss}.amps);
    end
  else
    tot=tot+sum(a.srccat.amps.*b.srccat.amps);
  end
  val=val+mpi_allreduce(tot);
end

