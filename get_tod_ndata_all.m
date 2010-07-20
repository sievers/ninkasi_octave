function[ndata]=get_tod_ndata_all(tods)
ndata=zeros(size(tods));
for j=1:length(tods),
  ndata(j)=get_tod_ndata(tods);
end

