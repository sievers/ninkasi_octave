function[ncut]=cut_tod_ends(tods,tcut)
ncut=zeros(length(tods),1);
for j=1:length(tods),
    ncut(j)=cut_tod_ends_c(tods(j),tcut);
end
