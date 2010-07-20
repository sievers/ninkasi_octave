function[mat]=get_darkdet_avg_cols(tod)

[darkdat,darkrow,darkcol]=get_dark_timestreams(tod);
darkcol=darkcol+1;
ncol=max(darkcol);

mat(size(darkdat,1),max(darkcol))=0;
for j=1:ncol,
  if sum(darkcol==j)>0,
    mat(:,j)=mean(darkdat(:,darkcol==j),2);
  end
end
