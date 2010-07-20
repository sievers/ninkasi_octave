function[map]=taper_map(map,template,fac)
if ~exist('fac')
  fac=0.1;
end

wtmax=max(max(template));
template=template/wtmax;
template=(fac-template)/(fac/10);
map=map./(1+exp(template));