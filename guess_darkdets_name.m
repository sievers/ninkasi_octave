function[darkname]=guess_darkdets_name(todname,dirroot)
darkname='';
if nargin<2
  dirroot='/home/sievers/act/darkDets/season2_cuts3/';
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
darkname=[dirroot '/' mydate '/' todname ext '/darkDets.txt'];
fid=fopen(darkname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end

darkname=[dirroot '/' mydate '/' todname ext '_darkDets.txt'];
fid=fopen(darkname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end



mydate=date_from_ctime(tok+86400,'mmm_dd_yyyy');
darkname=[dirroot '/' mydate '/' todname ext '/darkDets.txt'];
fid=fopen(darkname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end

darkname=[dirroot '/' mydate '/' todname ext '_darkDets.txt'];
fid=fopen(darkname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end


mydate=date_from_ctime(tok-86400,'mmm_dd_yyyy');
darkname=[dirroot '/' mydate '/' todname ext '/darkDets.txt'];
fid=fopen(darkname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end

darkname=[dirroot '/' mydate '/' todname ext '_darkDets.txt'];
fid=fopen(darkname,'r');
if (fid~=-1)
  fclose(fid);
  return;
end


darkname='';  %we didn't find anything here.
return
