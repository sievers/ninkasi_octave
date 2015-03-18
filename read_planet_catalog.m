function[cat]=read_planet_catalog(fname,ra_wrap)
lines=read_text_file(fname);
nl=length(lines);

for ii=nl:-1:1,
  ll=lines{ii};
  tags=strsplit(strtrim(ll),' ',true);
  cat.todname(ii)=tags(1);
  cat.uranus_ra(ii)=str2num(tags{8});
  cat.uranus_dec(ii)=str2num(tags{9});
end

if exist('ra_wrap')
  ind=cat.uranus_ra>ra_wrap;
  cat.uranus_ra(ind)=cat.uranus_ra(ind)-360;
end
