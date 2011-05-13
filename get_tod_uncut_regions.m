function[value]=get_tod_uncut_regions(tods)
for j=1:length(tods),
  get_tod_uncut_regions_c(tods(j));
  get_tod_cut_regions_c(tods(j));
  set_tod_cuts_global_indices_c(tods(j));
end
