function[ncut]=set_tod_window(tods,tcut)
ncut=zeros(length(tods),1);
for j=1:length(tods),
    ncut(j)=set_tod_window_c(tods(j),tcut);
end
