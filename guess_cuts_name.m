function[cutname]=guess_cuts_name(todname,dirroot)
cutname='';
if nargin<2
  %dirroot='/mnt/raid-cita/sievers/act/cuts/season2_cuts2/';
  %dirroot='/cita/d/raid-sievers3/sievers/act/cuts/season2_cuts2/';
  %dirroot='/home/sievers/act/cuts/season2_cuts2/';
  dirroot='/home/sievers/act/cuts/season2_cuts3/';
  %dirroot='/home/sievers/act/cuts/season2_cuts/';
end

if ~ischar(todname)
  assert(length(todname)==1);
  todname=get_tod_name(todname);
end
if todname(end)=='/'
  todname=todname(1:end-1);
end

%disp(todname)

[dr,todname,ext]=fileparts(todname);
tok=strtok(todname,'.');
tok=sscanf(tok,'%d');

mydate=date_from_ctime(tok,'mmm_dd_yyyy');
cutname=[dirroot '/' mydate '/' todname ext '/cuts.txt'];
fid=fopen(cutname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end

cutname=[dirroot '/' mydate '/' todname ext '_cuts.txt'];
fid=fopen(cutname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end


mydate=date_from_ctime(tok+86400,'mmm_dd_yyyy');
cutname=[dirroot '/' mydate '/' todname ext '/cuts.txt'];
fid=fopen(cutname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end



cutname=[dirroot '/' mydate '/' todname ext '_cuts.txt'];
fid=fopen(cutname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end

mydate=date_from_ctime(tok-86400,'mmm_dd_yyyy');
cutname=[dirroot '/' mydate '/' todname ext '/cuts.txt'];
fid=fopen(cutname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end


cutname=[dirroot '/' mydate '/' todname ext '_cuts.txt'];
fid=fopen(cutname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end


cutname='';  %we didn't find anything here.
return
