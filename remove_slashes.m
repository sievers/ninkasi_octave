function[tod_names]=remove_slashes(tod_names)
for j=1:length(tod_names),
  tt=tod_names{j};
  if tt(end)=='/', 
    tod_names(j)={tt(1:end-1)};
  end;
end
