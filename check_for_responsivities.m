function[have_responsivities]=check_for_responsivities(tod_names,dirroot)
have_responsivities=true(size(tod_names));
if ~exist('dirroot')
  dirroot='/home/sievers/act/responsivities/2009/';
end
for j=1:length(tod_names)
  fname=[dirroot '/' get_last_part(tod_names{j})];
  fid=fopen(fname,'r');
  if fid==-1
    have_responsivities(j)=false;
  else
    fclose(fid);
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[value]=get_last_part(tod_name)
ind=max(find(tod_name=='/'));
value=tod_name(ind+1:end);