function[offsets]=parse_tod_offsets(tod_offsets)
tod_ctimes=zeros(length(tod_offsets),1);
dx=tod_ctimes;
dy=tod_ctimes;
for j=1:length(tod_offsets),
  flub=strsplit(tod_offsets{j},' ',true);
  if strcmp(flub{1},'None')
    continue
  end

  fwee=strsplit(flub{1},'.',true);
  tod_ctimes(j)=str2num(fwee{1});
  dx(j)=str2num(flub{end-1});
  dy(j)=str2num(flub{end});
end
offsets.ctime=tod_ctimes;
offsets.dx=dx;
offsets.dy=dy;

