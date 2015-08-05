function[value]=cut_tod_azrange(tods,varargin)
if length(tods)>1
  for j=1:length(tods),
    cut_tod_azrange(tods(j),varargin{:});
  end
  return
end

assert(class(tods)=='int64');
azcut=get_keyval_default('azcut',[-62 -57;57 62; 72 76],varargin{:});
elthresh=get_keyval_default('elthresh',38,varargin{:});

[alt,az]=get_tod_altaz(tods);
alt=alt*180/pi;
az=az*180/pi;
az(az>180)=az(az>180)-360;
isbad=false(size(alt));
for j=1:size(azcut,1),
  ii=(alt<elthresh)&(az>=azcut(j,1))&(az<=azcut(j,2));
  isbad(ii)=true;
end

[istart,istop]=cutvec2inds(isbad);
for j=1:length(istart),
  cut_tod_global_c(tods,istart(j)-1,istop(j));
end


