function[x,b,d,Ad]=simple_pcg_octave_double(tods,x_in,b)

tmp=make_map_copy(x_in);
x=skymap2octave(x_in);


ntod=length(tods);
if ~exist('b')
  b=make_initial_rhs(tods,tmp);
end


tic;ax=mapset2mapset_octave(tods,x,tmp);toc;
r=b;
b=b-ax;
rr=sum(sum(r.^2));
r0sqr=rr;
iter=1;
tol=1e-8;
d=r;
disp(['r0sqr is ' num2str(r0sqr)]);
while (rr>r0sqr*tol)&(iter<20),
  Ad=mapset2mapset_octave(tods,d,tmp);
  dAd=sum(sum(d.*Ad));
  disp(['dAd is ' num2str(dAd)]);
  if (dAd==0)
    return;
  end
  alpha=rr/dAd;
  disp(['Alpha is ' num2str(alpha)]);
  x=x+d*alpha;
  r=r-Ad*alpha;
  rpp=sum(sum(r.^2));
  beta=rpp/rr;
  disp(['beta is ' num2str(beta)]);
  rr=rpp;
  d=r+beta*d;
  disp(num2str([iter rr]));
  iter=iter+1;
end


destroy_map(tmp);

return



function[b]=make_initial_rhs(tods,tmp)
ntod=length(tods);
clear-map(tmp);
for j=1:ntod,
  mytod=tods(j);
  allocate_tod_storage(mytod);
  read_tod_data(mytod);
  tic
  subtract_tod_median(mytod);
  toc
  tod2map(mytod,tmp);
  free_tod_storage(mytod);
end
b=skymap2octave(tmp);


function[cm]=subtract_tod_median(tod)
dat=get_tod_data(tod);
cm=omp_median(dat,2);
dat=dat-repmat(cm,[1 size(dat,2)]);
push_tod_data(dat,tod);
return
