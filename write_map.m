function[value]=write_map(map,fname,varargin)
if isstruct(map)
  if isfield(map,'partition')
    %write_partitioned_map(map,fname,varargin);
    fname=postpend_tag(fname,'.fits');
    write_partitioned_skymap_cea_fits(map,fname,varargin{:});
    return
  end
  
  octave2skymap(map);
  map=map.mapptr;
end
assert(class(map)=='int64');




mytype=get_map_type_c(map);
switch mytype
 case{'cea'}
  fname=postpend_tag(fname,'.fits');
  write_skymap_cea_fits(map,fname,varargin{:});
 case{'car'}
  fname=postpend_tag(fname,'.fits');
  write_skymap_car_fits(map,fname,varargin{:});
 case{'tan'}
  fname=postpend_tag(fname,'.fits');
  write_skymap_tan_fits(map,fname);
 case{'radec'}  
  %fname=postpend_tag(fname,'.map');
  %write_simple_map_c(map,fname);

  write_skymap_radec_fits(map,fname,varargin{:});
 case{'ring'}
  fname=postpend_tag(fname,'.fits');
  %mm=skymap2octave(map);
  %simple_write_healpix(fname,mm,[],'ring');
  write_healpix_ring(fname,map);
 case{'nest'}
  fname=postpend_tag(fname,'.fits');
  mm=skymap2octave(map);
  simple_write_healpix(fname,mm,[],'nest');
  
 otherwise
  disp('unrecognized type in write_map.');
end



function[fname]=postpend_tag(fname,tag)
%if there's no tag, add it to the end.
[a,b,c]=fileparts(fname);
if isempty(c)
  fname=[fname tag];
else
  if (length(c)>5)|(isdigit(c(1))),
    fname=[fname tag];
  end
end
