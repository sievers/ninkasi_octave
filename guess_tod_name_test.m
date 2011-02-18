function[nm_test]=guess_tod_name(todname,dirroot,force_find,strict)
if iscell(todname)
  if exist('dirroot')

    if ~exist('strict')
      strict=false;
    end
    if ~exist('force_find')
      force_find=false;
    end

    for j=length(todname):-1:1,
      nm_test{j}=guess_tod_name(todname{j},dirroot,force_find,strict);
    end
  else
    for j=length(todname):-1:1,
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
  %dirroot='/scratch/nolta/act/data/season2/merlin/'; %was this


  %dirroot='/project/nolta/act/data/season2/merlin/';
  dirroot=['/project/rbond/act/data/season' num2str(guess_tod_season(todname))  '/merlin/'];
  
end
if isempty(dirroot)
  %%dirroot='/scratch/sievers/act/data/season2/merlin/';
  %dirroot='/scratch/nolta/act/data/season2/merlin/';  %was this
  %dirroot='/project/nolta/act/data/season2/merlin/';
  dirroot=['/project/rbond/act/data/season' num2str(guess_tod_season(todname))  '/merlin/'];
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


nm_test=[dirroot mydate '/' tt ext];
if check_dirfile_exist(nm_test)
  return;
end


mydate=date_from_ctime(tok+86400,'yyyymmdd');


nm_test=[dirroot  mydate '/' tt  ext];
if check_dirfile_exist(nm_test)
  return;
end

nm_test=[dirroot  mydate '/' tt ext '/' tt ext ];
if check_dirfile_exist(nm_test)
  return;
end



nm_test='';
return






function[value]=check_strict(fname,tok)
disp(['checking ' fname '.zip'])
if exist([fname '.zip'],'file')==2
  value=true;
  return
end

if exist([fname '/Extras/' tok '_dat.run'],'file')
  value=true;
  return;
end
value=false;
return;   %remove above here to restore previous behavior

runname=[fname '.zip'];
fid=fopen(runname,'r');
if fid ~=-1
  value=true;
  fclose(fid);
  return;
end


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

function[value]=check_dirfile_exist(nm)
value=true;
if exist([nm '.zip'],'file')
  return;
end

if exist([nm '/format'],'file')
  return;
end

value=false;
return

