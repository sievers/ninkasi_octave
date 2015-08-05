function[i1,i2]=cutvec2inds(isbad)
mydiff=diff(isbad);
i1=find(mydiff>0)+1;
i2=find(mydiff<0);
if isbad(1)
  i1=[1;i1];
end
if isbad(end)
  i2=[i2;length(isbad)];
end

