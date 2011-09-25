function[mapset]=setup_mapset(tods,lims,myopts,outroot,projtype)

map=allocate_ninkasi_skymap(myopts.pixsize,lims(1)-myopts.pad,lims(2)+myopts.pad,lims(3)-myopts.pad,lims(4)+myopts.pad);
switch(projtype)
  case 'tan',
   set_map_projection_tan_simple(map);
 case 'cea',
  set_skymap_cea_simple_c(map);
 otherwise,
  error(['unsupported map projection type: ' projtype]);
  break;
end




mapset.skymap.mapptr=map;
mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);

if myopts.find_modes
  for j=1:length(tods),
    mdisp(['crunching  tod ' num2str(j)]);
    [rows,cols]=get_tod_rowcol(tods(j));
    mapset.corrnoise(j).map=zeros(get_tod_ndata(tods(j)),3);
    mapset.corrnoise(j).vecs=[ones(1,how_many_dets_are_kept(tods(j)));rows'-mean(rows);cols'-mean(cols)];
  end
end

