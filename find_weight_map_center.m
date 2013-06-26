function[xx,yy]=find_weight_map_center(map,varargin)
%find the center of the weight in a weightmap, for use in getting a preconditioner.  
%will find the centroid, then (optionally) search for the highest hit-count within a box around there in case it hit a zero.

boxsize=get_keyval_default('box_size',0,varargin{:});


xsum=cumsum(sum(map'));
xsum=xsum/xsum(end);
ysum=cumsum(sum(map));
ysum=ysum/ysum(end);


[a,xx]=min(abs(xsum-0.5));
[a,yy]=min(abs(ysum-0.5));

if boxsize>0,
  patch=map((xx-boxsize):(xx+boxsize),(yy-boxsize):(yy+boxsize));
  [a,dx]=max(patch);
  [a,dy]=max(a);dx=dx(dy);
  xx=xx-(boxsize+1)+dx;
  yy=yy-(boxsize+1)+dy;
  [dx dy a map(xx,yy)]
end



