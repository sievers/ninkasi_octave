function[found]=read_cuts(tods,dirroot)
found=true(size(tods));
for j=1:length(tods)
  if exist('dirroot')
    cutname=guess_cuts_name(tods(j),dirroot);
  else
    cutname=guess_cuts_name(tods(j));
  end
  if ~isempty(cutname),
    read_cuts_c(tods(j),cutname);
  else
    disp(['unable to find cuts on ' get_tod_name(tods(j))]);
    found(j)=false;
  end
  
end
