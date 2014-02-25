function[tods,tod_names]=cull_short_tods(tods,tod_names,len,dofree)
ind=true(size(tods));
if ~exist('dofree')
  dofree=false;
end

for j=1:length(tods),
  if get_tod_ndata(tods(j))<len
    ind(j)=false;
    if dofree
      erase_tod_c(tods(j));
    end

  end
end

tods=tods(ind);
if ~isempty(tod_names)
  tod_names=tod_names(ind);
end

if sum(~ind)
  disp(['culling ' num2str(sum(~ind)) ' short TODs']);
end
