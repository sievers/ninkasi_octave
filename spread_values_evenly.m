function[inds,sums]=spread_values_evenly(vals,n)
inds=zeros(size(vals));
sums(1:n,1)=0;
for j=1:length(vals),
  [a,b]=min(sums);
  sums(b)=sums(b)+vals(j);
  inds(j)=b;
end


