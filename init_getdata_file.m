function[myf,fmt]=init_getdata_file(todname)
if iscell(todname)
  for j=1:length(todname)
    myf(j)=init_getdata_file_c(todname{j});
    if nargout>1
      fmt(j)={read_dirfile_format(todname{j})};
    end
  end
else
  myf=init_getdata_file_c(todname);
  if nargout>1
    fmt=read_dirfile_format(todname);
  end
end

