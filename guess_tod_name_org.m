function[nm_test]=guess_tod_name(todname,dirroot,force_find,strict)
if iscell(todname)
  if exist('dirroot')

    if ~exist('strict')
      strict=false;
    end
    if ~exist('force_find')
      force_find=false;
    end

    for j=1:length(todname),
      nm_test{j}=guess_tod_name(todname{j},dirroot,force_find,strict);
    end
  else
    for j=1:length(todname),
      nm_test{j}=guess_tod_name(todname{j});
    end
  end
  
  return
end
if ~exist('force_find')
  force_find=false;
end

if ~exist('strict')
  strict=false;
end


if nargin<2
  %dirroot='/mnt/raid-cita/sievers/act/cuts/season2_cuts2/';
  %dirroot='/cita/d/raid-sievers3/sievers/act/cuts/season2_cuts2/';
  %dirroot='/scratch/sievers/act/data/season2/merlin/';
  dirroot='/scratch/nolta/act/data/season2/merlin/';
end
if isempty(dirroot)
  %dirroot='/scratch/sievers/act/data/season2/merlin/';
  dirroot='/scratch/nolta/act/data/season2/merlin/';
end

[dr,tt,ext]=fileparts(todname);

tok=strtok(tt,'.');

tok_org=tok;
tok=sscanf(tok,'%d');

mydate=date_from_ctime(tok,'yyyymmdd');
nm_test=[dirroot  mydate '/' tt ext '/'];
%disp(['starting nm_test is ' nm_test]);
if force_find
  return;
end



fid=fopen([nm_test 'format'],'r');
if (fid~=-1)
  fclose(fid);
  if strict
    if ~check_strict(nm_test,tok_org)
      nm_test='';
    end
  end
  return;
end


mydate=date_from_ctime(tok+86400,'yyyymmdd');
nm_test=[dirroot  mydate '/' tt  ext '/'];
fid=fopen([nm_test 'format'],'r');
if (fid~=-1)
  fclose(fid);
  if strict
    if ~check_strict(nm_test,tok_org)
      nm_test='';
    end
  end
  return;
end



nm_test=[dirroot  mydate '/' tt ext '/' tt ext '/'];
fid=fopen([nm_test 'format'],'r');
if (fid~=-1)
  fclose(fid);
  if strict
    if ~check_strict(nm_test,tok_org)
      nm_test='';
    end
  end  
  return;
end




nm_test=[dirroot  mydate '/' tt ext '/' tt ext '/'];
fid=fopen([nm_test 'format'],'r');

if (fid~=-1)
  fclose(fid);
  if strict
    if ~check_strict(nm_test,tok_org)
      nm_test='';
    end
  end
  
  return;
end


 nm_test='';
return






function[value]=check_strict(fname,tok)
runname=[fname '/Extras/' tok '_dat.run'];
fid=fopen(runname,'r');
if fid==-1
  value=false;
  disp(['failed extra check on ' runname]);
else
  %disp(['passed extra check on ' runname]);
  value=true;
  fclose(fid);
end
