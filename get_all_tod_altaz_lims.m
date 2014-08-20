function[lims,mylims]=get_all_tod_altaz_lims(tods)
mylims=zeros(length(tods),4);
for j=1:length(tods),
  mylims(j,:)=get_tod_altaz_lims_c(tods(j));
end
mylims=[min(mylims(:,1)) max(mylims(:,2)) min(mylims(:,3)) max(mylims(:,4))];
lims=mpi_reduce_limits(mylims);
