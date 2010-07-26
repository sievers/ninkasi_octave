function[cuts]=find_cuts_from_weights(wt)
cut=9.9e5;
crud=diff(wt>cut);
ind=find(crud~=0);
ind=[0;ind;length(wt)];
if min(wt(ind(1)+1:ind(2)))>cut,
  jstart=1;
else
  jstart=2;
end
jj=1;
cuts=zeros(ceil(length(ind)/2),2);
for j=jstart:2:length(ind)-1,
  cuts(jj,1)=ind(j)+1;
  cuts(jj,2)=ind(j+1);
  jj=jj+1;
end
if size(cuts,1)==1
  if cuts(1,1)==1
    if cuts(1,2)==length(wt)
      cuts(1,2)=2^31;
      %disp('cutting now.');
    end
  end
end


cuts=cuts-1;



