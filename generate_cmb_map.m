function[map,rr]=generate_cmb_map(mapsize,pixsize,powspec)
if ~exist('powspec')
    powspec=wmap_5year_spec;
    powspec=powspec(:,2);
end
l=(1:length(powspec))';
powspec=powspec*2*pi./(l.*(l+1));



nx=get_smallest_factor_len(2*mapsize(1)+100,100);
ny=get_smallest_factor_len(2*mapsize(2)+100,100);



map=randn(nx,ny);
map=fft2(map);
rx=get_vec(nx);
ry=get_vec(ny);
pixsize_x=1/(pixsize/60*pi/180*nx)*2*pi;
pixsize_y=1/(pixsize/60*pi/180*ny)*2*pi;

rr=sqrt(repmat(pixsize_x^2*rx'.^2,[1 ny])+repmat(pixsize_y^2*ry.^2,[nx 1]));

%rr=interp1(l,sqrt(powspec),rr,'pchip',0);
rr=interp1(l,sqrt(powspec),rr,0);
map=ifft2(map.*rr)*sqrt(numel(rr))/2;

map=map(1:mapsize(1),1:mapsize(2));

%pixsize in arcmin



function[vec]=get_vec(nx);
if (iseven(nx)),
    vec=(-nx/2:(nx/2-1));
else, 
    vec=((-nx+1)/2):((nx-1)/2);
end;
vec=ifftshift(vec);