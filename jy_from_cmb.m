function[jy,domega]=jy_from_cmb(nu,arcsec)
if ~exist('nu')
    nu='';
end
if isempty(nu)
    nu=31;
end

k= 1.3806503e-16;
h=6.626068e-27;
T=2.725;
c=2.99792458e10;

x=h*(nu*1e9)/k/T;
fac=2*k/c^2*(k*T/h)^2*1e23;
amp=fac*x^4*exp(x)/(exp(x)-1)^2;
if exist('arcsec')
    domega=(arcsec/3600*pi/180)^2;
else
    domega=1;  %1 steradian;
end

jy=domega*amp;

