function[offsets]=parse_tod_offsets(tod_offsets)
tod_ctimes=zeros(length(tod_offsets),1);
dx=tod_ctimes;
dy=tod_ctimes;
arr=tod_ctimes;
for j=1:length(tod_offsets),
  flub=strsplit(strtrim(tod_offsets{j}),' ',true);
  if strcmp(flub{1},'None')
    continue
  end

  fwee=strsplit(strtrim(flub{1}),'.',true);
  tod_ctimes(j)=str2num(fwee{1});
  dx(j)=str2num(flub{end-1});
  dy(j)=str2num(flub{end});
  ff=fwee{end};
  arr(j)=str2num(ff(end));
end
offsets.ctime=tod_ctimes;
offsets.dx=dx;
offsets.dy=dy;
offsets.arr=arr;
