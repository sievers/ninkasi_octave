function[vec]=get_sawtooth(speed,amp,nsamp)
vec=(1:nsamp)*speed*pi/2/amp;
vec=2*amp*(acos(cos(vec'))/pi-0.5);

