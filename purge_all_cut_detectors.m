function[value]=purge_all_cut_detectors(tods)
for j=1:length(tods),
  purge_cut_detectors(tods(j));
end

