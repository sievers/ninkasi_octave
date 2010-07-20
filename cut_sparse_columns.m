function[ncut]=cut_sparse_columns(tods,nmin)
if ~exist('nmin')
  nmin=4;
end

if length(tods)>1,
  for j=1:length(tods),
    ncut(j)=cut_sparse_columns(tods(j),nmin);
  end
  return;
end
[rows,cols]=get_tod_rowcol(tods);
cc=unique(cols);
ncut=0;
to_cut=[];
for j=1:length(cc),
  if sum(cols==cc(j))<nmin
    to_cut=[to_cut cc(j)];
    ncut=ncut+sum(cols==cc(j));
  end
end
cut_column(tods,to_cut);
      

