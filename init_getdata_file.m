function[myf]=init_getdata_file(todname)
if iscell(todname)
  for j=1:length(todname)
    myf(j)=init_getdata_file_c(todname{j});
  end
else
  myf=init_getdata_file_c(todname);
end

