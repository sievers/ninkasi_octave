function[params,wt]=fit_loglog_poly(y,x,order,wt_scheme);
if size(x,1)==1,
  x=x';
end
if size(y,1)==1,
  y=y';
end

df=diff(x);
wt1=df./x(1:end-1);
wt2=df./x(2:end);
wt=0.5*([wt1(1);wt1] + [wt2;wt2(end)]);

x(end+1)=x(end);
ii=round(0.75*length(y));
y(end+1)=median(y(ii:end));  %we're going to peg the end
disp(y(end));
wt(end+1)=sum(wt);

xmat=zeros(length(x),order);
xmat(:,1)=1;
logx=log(x);
for j=2:order,
  xmat(:,j)=xmat(:,j-1).*logx;
end

params=linfit(log(y),xmat,1./wt);
