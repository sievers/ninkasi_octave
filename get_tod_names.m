function[tod_names]=get_tod_names(tods)

tod_names={};


for j=1:length(tods),
  tod_names(j)={get_tod_name(tods(j))};
end

