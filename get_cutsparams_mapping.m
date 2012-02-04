function[vec]=get_cutsparams_mapping(tod,dt_per_param)
vec=(0:get_tod_ndata(tod)-1)*get_tod_dt(tod);
vec=ceil(vec/dt_per_param);

vec(1)=1;
vec(2:3)=2;
vec(4:6)=3;
vec(7:20)=4;
if length(vec)>20,
  vec(21:end)=vec(21:end)+4;
end
