function[ra,dec]=get_radec_from_loic_files(fname)
[a,b]=system(['awk ''{if ($6==2) print $0}'' < ' fname ]);
tod_names=strsplit(strtrim(b),sprintf('\n'),true);
ra=zeros(size(tod_names));
dec=0*ra;

for j=1:length(tod_names),
  tags=strsplit(strtrim(tod_names{j}),[' ' sprintf('\t')],true);
  tt=strsplit(strtrim(tags{1}),'.',true);
  ct=str2num(tt{1});
  alt=str2num(tags{3})*pi/180;
  az=str2num(tags{4})*pi/180;
  [myra,mydec]=get_radec_from_altaz_actpol_c(az,alt,ct);
  ra(j)=myra;
  dec(j)=mydec;
end

  
