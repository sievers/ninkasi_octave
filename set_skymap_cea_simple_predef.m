function[value]=set_skymap_cea_simple_predef(map,pixsize_in,pv)
if ischar(pixsize_in)
  found=0;
  if strcmp(pixsize_in,'act_equ')
    pixsize=0.00825;
    pv=1.0;
    found=1;
  end
  if strcmp(pixsize_in,'act_south')
    pixsize=0.013869678;
    pv=0.35381412;
    found=1;
  end
  assert(found==1);  %if this fails, we didn't have a valid tag for the settings
else
  pixsize=pixsize_in;
end
set_skymap_cea_simple_predef_c(map,pixsize,pv);
    
