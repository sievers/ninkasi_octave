function[dat]=recalibrate_data_pairs(dat,match_list)
for j=1:length(match_list),
  fwee=[dat(:,j) dat(:,match_list(j))];
  [vv,ee]=eig(fwee'*fwee);
  ee=diag(ee);
  [a,b]=max(ee);
  facs=vv(:,b);facs=facs/mean(facs);
  dat(:,j)=dat(:,j)/facs(1);
  dat(:,match_list(j))=dat(:,match_list(j))/facs(2);
end

  
