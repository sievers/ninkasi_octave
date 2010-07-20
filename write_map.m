function[value]=write_map(map,fname)
if isstruct(map)
  octave2skymap(map);
  map=map.mapptr;
end
assert(class(map)=='int64');




mytype=get_map_type_c(map);
switch mytype
 case{'cea'}
  fname=postpend_tag(fname,'.fits');
  write_skymap_cea_fits(map,fname);
 case{'radec'}  
  fname=postpend_tag(fname,'.map');
  write_simple_map_c(map,fname);
 otherwise
  disp('unrecognized type in write_map.');
end



function[fname]=postpend_tag(fname,tag)
%if there's no tag, add it to the end.
[a,b,c]=fileparts(fname);
if isempty(c)
  fname=[fname tag];
else
  if (length(c)>5)|(isdigit(c(1)))
    fname=[fname tag];
  end
end