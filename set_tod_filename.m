function[value]=set_tod_filename(tod,todname)
if iscell(todname)
  big_name=todname{1};
  for j=2:length(todname)
    big_name=sprintf('%s\n%s',big_name,todname{j});
  end  
  set_tod_filename_c(tod,big_name);
else
  set_tod_filename_c(tod,todname);
end
