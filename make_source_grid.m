function[map]=make_source_grid(n,m,dx,dy,width)

map=zeros(n,m);
xx=round(dx/2):dx:n;
yy=round(dy/2):dy:m;
for j=1:length(xx),
  map(xx(j),yy)=1;
end
map=smooth_image(map,width);
