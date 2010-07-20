function[value]=reverse_tod_uncut_regions(tods)
assert(class(tods)=='int64');
for j=1:length(tods),
  reverse_tod_uncut_regions_c(tods(j));
end
