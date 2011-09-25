function[mapset]=read_mapset(tods,dirname)
myid=mpi_comm_rank+1;
dirname(end+1)='/';
if (myid==1)
  try
    map=fits_image_read([dirname 'skymap.fits']);
  catch
    map=[];
  end
  mpi_bcast_array(map);
else
  map=mpi_bcast_array();
end
if ~isempty(map)
  mapset.skymap.map=map;
end
for j=1:numel(tods),
  tod_name=get_tod_names(tods(j));
  tag=get_tod_tags_from_names(tod_name);

  %do cutvecs
  fid=fopen([dirname 'cutvecs_' tag]);
  if (fid~= -1)
    n=fread(fid,1,'int');
    vec=fread(fid,[n 1],'double');
    fclose(fid);
    if numel(tods)==1
      mapset.cutvecs=vec;
    else
      mapset.cutvecs(j)={vec};
    end
  end
end

    