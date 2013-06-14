function[value]=save_mapset(mapset,tods,dirname)
myid=mpi_comm_rank+1;
if (myid==1)
  system(['mkdir ' dirname ' >& /dev/null']);
end
dirname(end+1)='/';
mpi_barrier;  %make sure everyone waits for mkdir to succees
if (myid==1)
  if isfield(mapset,'skymap')
    octave2skymap(mapset.skymap);
    write_map(mapset.skymap.mapptr,[dirname 'skymap']);
  end
end
if isfield(mapset,'cutvecs')  
  for j=1:length(tods),
    tod_name=get_tod_name(tods(j));
    if iscell(tod_name)
      tod_name=tod_name{1};
    end
    tag=get_tod_tags_from_names(tod_name);
    fid=fopen([dirname 'cutvecs_' tag],'w');
    if iscell(mapset.cutvecs)
      fwrite(fid,numel(mapset.cutvecs{j}),'int');
      fwrite(fid,mapset.cutvecs{j},'double');
    else
      fwrite(fid,numel(mapset.cutvecs),'int');
      fwrite(fid,mapset.cutvecs,'double');
    end

    fclose(fid);
  end
end


  
