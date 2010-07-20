function[ndet]=get_tod_ndet(tods)
ndet=zeros(size(tods));
for j=1:length(tods),
  [rows,cols]=get_tod_rowcol(tods(j));
  ndet(j)=numel(rows);
end
