function[map]=center_image(map)
[a,b]=max(abs(map));
[a,c]=max(a);
b=b(c);
%[a map(b,c)]

if (b>ceil(size(map,1)/2)),
  map(2*b,1)=0;
else
  npad=size(map,1)-2*b;
  map=[zeros(npad,size(map,2));map];
end



if (c>ceil(size(map,2)/2)),
  map(1,2*c)=0;
else
  npad=size(map,2)-2*c;
  map=[zeros(size(map,1),npad) map];
end

map=insert_map_into_map(zeros(max(size(map))),map);

  
