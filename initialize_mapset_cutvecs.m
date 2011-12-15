function[cutvecs]=initialize_mapset_cutvecs(tods)
cutvecs=cell(size(tods));
for j=1:length(tods),
  cutvecs(j)={zeros(get_numel_cut_c(tods(j)),1)};
end

