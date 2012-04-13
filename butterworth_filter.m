function[filt]=butterworth_filter(n,dt,expt)

if ~exist(expt)
  expt='act';
end

if isequal(expt,'act')
  b_1_1=(-2.0*32092.0/2^15);
  b_1_2=( 2.0*15750.0/2^15);
  b_2_1=(-2.0*31238.0/2^15);
  b_2_2=( 2.0*14895.0/2^15);
  f_samp=1./(0.00000002*100.*33.);
  mynorm=4/(1+b_1_1+b_1_2)*4/(1+b_2_1+b_2_2);
  omega=(0:n-1)';
  omega=omega/dt/n/f_samp*2*pi;
  h1_omega=(1+2*exp(-i*omega)+exp(-2*i*omega))./(1+b_1_1*exp(-i*omega)+b_1_2*exp(-2*i*omega));
  h2_omega=(1+2*exp(-i*omega)+exp(-2*i*omega))./(1+b_2_1*exp(-i*omega)+b_2_2*exp(-2*i*omega));
  filt=h1_omega.*h2_omega/mynorm;
  return
end
if isequal(expt,'abs')
  b_1_1=-2*32022/2^15;
  b_1_2=2*15648/2^15;
  b_2_1=-2*32449/2^15;
  b_2_2=2*16075/2^15;
  
  f_samp=1./(0.00000002*150.*22.);
  mynorm=4/(1+b_1_1+b_1_2)*4/(1+b_2_1+b_2_2);
  omega=(0:n-1)';
  omega=omega/dt/n/f_samp*2*pi;
  h1_omega=(1+2*exp(-i*omega)+exp(-2*i*omega))./(1+b_1_1*exp(-i*omega)+b_1_2*exp(-2*i*omega));
  h2_omega=(1+2*exp(-i*omega)+exp(-2*i*omega))./(1+b_2_1*exp(-i*omega)+b_2_2*exp(-2*i*omega));
  filt=h1_omega.*h2_omega/mynorm;


end
  

