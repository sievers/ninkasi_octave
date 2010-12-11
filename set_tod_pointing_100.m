function[lims,isrising]=set_tod_pointing_100(tod,my_poff)

my_delts=get_tod_pointing_offsets_type_100(tod,my_poff);
%mdisp('got offsets');
my_dalt=-111*ones(my_poff.nrow,my_poff.ncol);
my_daz=-111*ones(my_poff.nrow,my_poff.ncol);

for k=1:length(my_delts),
  my_daz(my_delts(k,1)+1,my_delts(k,2)+1)=my_delts(k,3);
  my_dalt(my_delts(k,1)+1,my_delts(k,2)+1)=my_delts(k,4);
end

%mdisp('setting pointing c');
lims=set_tod_pointing_c(tod,my_dalt,my_daz);
%mdisp('getting altaz');
[myaz,myalt,myctime]=get_median_altaz_c(tod);
if myaz<pi
  isrising=true;
else
  isrising=false;
end
