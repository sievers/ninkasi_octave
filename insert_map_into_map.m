function[map,mybox1,mybox2]=insert_map_into_map(map,map_in)
map(:,:)=0;

s1=size(map);
s2=size(map_in);
[xa1,xa2,xb1,xb2]=inbounds_sizes(s1(1),s2(1));
[ya1,ya2,yb1,yb2]=inbounds_sizes(s1(2),s2(2));

mybox1=[xa1 xa2 ya1 ya2];
mybox2=[xb1 xb2 yb1 yb2];
map(xa1:xa2,ya1:ya2)=map_in(xb1:xb2,yb1:yb2);


function[xa1,xa2,xb1,xb2]=inbounds_sizes(n1,n2)
if n2>n1
  xa1=1;
  xa2=n1;

  xb1=1+(n2-n1);
  xb2=n2;
else
  dx=round((n1-n2)/2);
  xa1=dx+1;
  xa2=dx+n2;
  xb1=1;
  xb2=n2;
end
