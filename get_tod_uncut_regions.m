function[value]=get_tod_uncut_regions(tods)
for j=1:length(tods),
  get_tod_uncut_regions_c(tods(j));
end
