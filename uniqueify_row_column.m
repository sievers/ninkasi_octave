function[rr2,cc2]=uniqueify_row_column(rr,cc)

vec=rr+i*(cc+1);
if length(unique(vec))==numel(vec), return;end; %nothing to do here

cc2=cc;rr2=rr;


cols=unique(cc);
for mycol=1:length(cols),
  ind=(cc==cols(mycol));
  tmp=rr(ind);
  if length(unique(tmp))<numel(tmp),
    tmp=(1:length(tmp))';
    rr2(ind)=tmp;
  end
end
