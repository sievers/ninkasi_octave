function[value]=set_skymap_car_simple_predef(map,decpix,rapix)
if ischar(decpix)
  found=0;
  if strcmp(decpix,'act_equ')
    decpix=0.00825;
    rapix=-decpix;
    found=1;
  end
  if strcmp(pixsize_in,'act_south')
    decpix=0.013869678;
    pv=0.35381412;
    rapix=-decpix/sqrt(px); %roughly  correct I think
    found=1;
  end
  assert(found==1);  %if this fails, we didn't have a valid tag for the settings
else
  if ~(exist('rapix')==1),
    rapix=-decpix;
  end
end
set_skymap_car_simple_predef_c(map,decpix,rapix);
    
