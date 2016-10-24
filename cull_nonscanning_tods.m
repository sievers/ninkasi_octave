function[tods,tod_names]=cull_nonscanning_tods(tods,tod_names,thresh,dofree)
if ~exist('thresh')
  thresh=1e-4;  %some sensible value for the average diff of the az
end

ind=true(size(tods));
if ~exist('dofree')
  dofree=true;
end

for j=1:length(tods),
  %if get_tod_ndata(tods(j))<len
  [alt,az]=get_tod_altaz(tods(j));
  if (mean(abs(diff(unwrap(az))))<thresh),
    ind(j)=false;
    disp(['culling non-scanning tod ' tod_names{j}])
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
