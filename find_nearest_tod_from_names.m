function[match,dt]=find_nearest_tod_from_names(nm,nm2)
t1=get_tod_ctimes_from_names(nm);
t2=get_tod_ctimes_from_names(nm2);
match=cell(size(nm));
dt=zeros(length(nm),1);
for j=1:length(nm),
  [a,b]=min(abs(t2-t1(j)));
  match(j)=nm2(b);
  dt(j)=a;
end