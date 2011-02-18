function[det_offs]=get_tod_pointing_offsets_type_100(tod,poff)
[pp,myalt]=select_pointing_block(tod,poff.pstruct);

offsets=calculate_boresight_offsets(tod,pp);
det_offs=poff.det_raw;
det_offs(:,3)=det_offs(:,3)+offsets(1);
det_offs(:,4)=det_offs(:,4)+offsets(2);
det_offs(:,5:6)=det_offs(:,3:4);  %stick in x/y offsets as well
[dalt,daz]=convert_offsets_to_nomdelts(det_offs(:,3),det_offs(:,4),myalt);

det_offs(:,3)=daz;
det_offs(:,4)=dalt;

return



function[offsets]=calculate_boresight_offsets(tod,pp)
offsets=zeros(2,1);
if isfield(pp,'apply_shift')
  offsets=offsets+pp.apply_shift;
end
if isfield(pp,'apply_dy_ctime')
  myctime=get_tod_ctime_c(tod);
  dt=myctime-pp.apply_dy_ctime(1);
  flipud(pp.apply_dy_ctime(2:end));
  offsets(2)=offsets(2)+polyval(flipud(pp.apply_dy_ctime(2:end)),dt);
end

%do a time-of-day correction
if isfield(pp,'apply_dy_dtime')
  myctime=get_tod_ctime_c(tod);
  dt=rem((myctime-pp.apply_dy_dtime(1)),86400);
  flipud(pp.apply_dy_ctime(2:end));
  delta=polyval(flipud(pp.apply_dy_dtime(2:end)),dt);
  mdisp(['dtime delta on ' get_tod_name(tod) ' is ' num2str(delta)]);
  offsets(2)=offsets(2)+delta;
end


if isfield(pp,'apply_dx_ctime')
  myctime=get_tod_ctime_c(tod);
  dt=myctime-pp.apply_dx_ctime(1);
  offsets(1)=offsets(1)+polyval(flipud(pp.apply_dx_ctime(2:end)),dt);
end



function[pp,myalt]=select_pointing_block(tod,pstructs)
[myaz,myalt,myctime]=get_median_altaz_c(tod);

for j=1:length(pstructs)
  pp=pstructs{j};
  do_i_match=match_ctime(myctime,pp);
  do_i_match=do_i_match&(match_az(myaz,pp));

  if (do_i_match)
    return;
  end
end
%pp=[];  %if here, it's 'cause we didn't match any of the candidate blocks
error(['matching pointing block not found on ' get_tod_name(tod) ' with ' num2str([myctime myaz])]);
return

function[tf]=match_ctime(myctime,pp)
tf=true;
if isfield(pp,'select_ctime') %if no ctime specified, default to true
  if myctime<pp.select_ctime(1)
    tf=false;
  end
  if myctime>pp.select_ctime(2)
    tf=false;
  end
end
return

function[tf]=match_az(myaz,pp)
tf=true;
if isfield(pp,'select_az')
  if myaz*180/pi<pp.select_az(1)
    tf=false;
  end
  if myaz*180/pi>pp.select_az(2)
    tf=false;
  end
end
return