function[alt,az]=get_altaz_from_radec_ctime(ra,dec,ct)
if numel(ct)==1
  ct=0*ra+ct;
end


assert(numel(ra)==numel(dec));

if numel(ra)==1
  ra=repmat(ra,size(ct));
  dec=repmat(dec,size(ct));
end
assert(numel(ra)==numel(ct));

alt=0*ra;
az=0*ra;


for j=1:length(ra)
  [alt(j),az(j)]=get_altaz_from_radec_ctime_c(ra(j),dec(j),ct(j));
  %alt(j)=tmpalt;
  %az(j)=tmpaz;
end
