function[value]=photometer_point(map,x,y,r1,r2)
%return the total flux in a circle centered on x and y, within radius r1
%after subtracting off the median in the donut between r1 and r2
if ~exist('r1')
  r1=5;
end
if ~exist('r2')
  r2=15;
end


x=round(x);
y=round(y);
rr=ceil(r2);
mm=map([x-rr:x+rr],[y-rr:y+rr]);
vec=-rr:rr;
rmat=bsxfun(@hypot,vec,vec');
m1=mm(rmat<=r1);
m2=mm((rmat<=r2)&(rmat>r1));
value=sum(m1)-median(m2)*numel(m1);


% 03:51:31.066  -51:43:44.918 x=4424; y=962;      1.2822e+04   1.6723e+04   1.2259e+04   1.4613e+04   1.4277e-01   3.1648e-01
% 03:14:28.361  -51:04:10.700 x=5084; y=1042;     9.8391e+03   1.0629e+04   9.7317e+03   1.1199e+04   6.6981e-02   1.4176e-01
% 04:25:08.650  -53:32:33.461 x=3823; y=739;      1.8294e+04   1.8554e+04   1.9276e+04   1.7446e+04   4.1079e-02   9.9507e-02
% 05:40:48.214  -54:18:31.486 x=2473; y=648;      4.4882e+04   4.3706e+04   4.4593e+04   4.4685e+04   1.1721e-02   2.6448e-02
% 05:49:46.236  -52:46:54.966 x=2313; y=832;      1.6591e+04   1.6837e+04   1.7250e+04   1.5923e+04   3.3391e-02   7.9711e-02
