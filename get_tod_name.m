function[nm]=get_tod_name(tod)
nm=get_tod_name_c(tod);
cc=sprintf('\n');
nm=strsplit(strtrim(nm),cc,true);
if numel(nm)==1
  nm=nm{1};
end

