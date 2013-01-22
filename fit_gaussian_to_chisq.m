function[value]=fit_gaussian_to_chisq(map,xvec,yvec)
if ~exist('xvec')
  xvec=(1:size(map,1))';
end
if ~exist('yvec')
  yvec=xvec;
end

[a,b]=min(map);
[a,c]=min(a);
b=b(c);
pad=2;

i1=b-pad;
i2=b+pad;

if i1<1
  i1=1;
end
if i2>size(map,1)
  i2=size(map,1);
end

fitp=polyfit(xvec(i1:i2)',map(i1:i2,c),2);
xcent=-0.5*fitp(2)/fitp(1);

i1=c-pad;
i2=c+pad;

if i1<1
  i1=1;
end
if i2>size(map,2)
  i2=size(map,2);
end

fitp=polyfit(yvec(i1:i2),map(b,i1:i2),2);
ycent=-0.5*fitp(2)/fitp(1);
value=[xcent ycent];

