function[value]=find_happy_fft_pad(val,pvec)
if ~exist('pvec')
  pvec=[2 3 5];
end
if numel(val)>1
  for j=1:length(val),
    value(j)=find_happy_fft_pad(val(j),pvec);
  end
  return
end

facs=[];
%first pull out any exact factors we have
for j=1:length(pvec),
  while rem(val,pvec(j))==0
    facs=[facs pvec(j)];
    val=val/pvec(j);
  end
end
vv=val;
[val vv]

%now just brute-force things
while max(factor(vv))>max(pvec),
  vv=vv+1;
end
value=vv;
value=value*prod(facs);

