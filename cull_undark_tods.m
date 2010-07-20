function[tods,tod_names]=cull_undark_tods(tods,tod_names,ndet_min,dofree)
ind=true(size(tods));
if ~exist('dofree')
  dofree=false;
end

for j=1:length(tods),
  [row,col]=get_darkdets_rowcol(tods(j));
  if length(row)<ndet_min
    ind(j)=false;
    if dofree
      erase_tod_c(tods(j));
    end
    
  end
end

tods=tods(ind);
tod_names=tod_names(ind);
if sum(~ind)
  disp(['culling ' num2str(sum(~ind)) ' undark TODs']);
end
