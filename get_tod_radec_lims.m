function[lims,ll]=get_tod_radec_lims(tods)
for j=length(tods):-1:1,
  ll(j,:)=get_tod_radec_lims_c(tods(j));
end
lims=[min(ll(:,1)) max(ll(:,2)) min(ll(:,3)) max(ll(:,4))];
