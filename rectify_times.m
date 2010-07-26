function[tt]=rectify_times(tt)

dd=diff(tt);
rat=max(diff(tt))/min(diff(tt));
if (rat < 1.05 ) & (rat>0.95)   %at this point, times are OK. seemingly backwards is if negative dt's happen.
  return;
end


tsamp=median(dd);
thresh=20;  %correct any samples where time is off by more than 5 times median scatter
scat=median(abs(dd-tsamp));
ind=abs(dd-tsamp)>thresh*scat;
%disp(['flagging ' num2str(sum(ind)) ' samples.'])
ind=find(ind);
extras=[];
for j=1:length(ind)-1,
  dind=ind(j+1)-ind(j);
  if (dind>1)&(dind<10) %we're missing samples
    extras=[extras (ind(j)+1):(ind(j+1)-1)];
  end
end
ind=sort([ind; extras']);
i2=find(diff(ind)>1);
i2=[1 ;i2; length(ind)];

for j=1:length(i2)-1,
  seg=ind(i2(j)+1:i2(j+1));
  tt(seg)=tt(seg(1)-1)+(1:length(seg))/(length(seg)+1)*(tt(seg(end)+1)-tt(seg(1)-1));
end

%for j=1:length(ind)
%  tt(ind(j))=tt(ind(j)-1)+tsamp;
%end

