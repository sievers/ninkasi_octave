function[b,d,Ad]=simple_pcg_octave(tods,x)

b=make_map_copy(x);

clear_map(b);
clear_map(x);

ntod=length(tods);
for j=1:ntod,
  mytod=tods(j);
  allocate_tod_storage(mytod);
  tic
  read_tod_data(mytod);
  toc
  tod2map(mytod,b);
  free_tod_storage(mytod);
end

tic;ax=mapset2mapset(tods,x);toc;
r=make_map_copy(b);
map_axpy(b,ax,-1);
rr=map_times_map(r,r);
r0sqr=rr;
iter=1;
tol=1e-8;
d=make_map_copy(r);
disp(['r0sqr is ' num2str(r0sqr)]);
while (rr>r0sqr*tol)&(iter<20),
  Ad=mapset2mapset(tods,d);
  dAd=map_times_map(d,Ad);
  disp(['dAd is ' num2str(dAd)]);
  if (dAd==0)
    return;
  end
  alpha=rr/dAd;
  disp(['Alpha is ' num2str(alpha)]);
  map_axpy(x,d,alpha);
  map_axpy(r,Ad,-alpha);
  rpp=map_times_map(r,r);
  beta=rpp/rr;
  disp(['beta is ' num2str(beta)]);
  rr=rpp;
  dd=make_map_copy(r);
  map_axpy(dd,d,beta);
  destroy_map(d);
  d=dd;
  destroy_map(Ad);
  disp(num2str([iter rr]));
  iter=iter+1;
end
destroy_map(b);
destroy_map(r);
destroy_map(d);
if (nargout<1)
  destroy_map(rhs);
end



  


%function[map2]=mapset2mapset(tods,map)
%map2=make_map_copy(map);
%clear_map(map2);
%for j=1:length(tods),
%  mytod=tods(j);
%  allocate_tod_storage(mytod);
%  map2tod(map,mytod);
%  tod2map(mytod,map2);
%end

