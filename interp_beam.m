function[new_beam]=interp_beam(beam,fac)
if ~exist('fac')
  fac=4;
end
dx=median(diff(beam(:,1)));
x=(0:(fac*length(beam)))'*dx/fac;
y=interp1(beam(:,1),beam(:,2),x,'spline');
ii=isfinite(y);
new_beam=[x(ii) y(ii)];
