function[x,b,weights]=simple_pcg_octave_double_precon(tods,x_in,b)

tmp=make_map_copy(x_in);
x=skymap2octave(x_in);

make_weightmap_octave(tods,tmp);
weights=skymap2octave(tmp);
clear_map(tmp);


ntod=length(tods);
if ~exist('b')
  b=make_initial_rhs(tods,tmp);
end


tic;ax=mapset2mapset_octave(tods,x,tmp);toc;
r=b-ax;
Mr=apply_precon(r,weights);
%b=b-ax;

iter=1;
tol=1e-8;
d=Mr;
rMr=sum(sum(r.*d));
r0sqr=rMr;
disp(['r0sqr is ' num2str(r0sqr)]);
while (rMr>r0sqr*tol)&(iter<20),
  Ad=mapset2mapset_octave(tods,d,tmp);
  dAd=sum(sum(d.*Ad));
  %disp(['dAd is ' num2str(dAd)]);
  if (dAd==0)
    return;
  end
  alpha=rMr/dAd;
  %disp(['Alpha is ' num2str(alpha)]);
  x=x+d*alpha;
  r=r-Ad*alpha;
  Mr=apply_precon(r,weights);
  rpp=sum(sum(r.*Mr));  
  beta=rpp/rMr;
  %disp(['beta is ' num2str(beta)]);
  rMr=rpp;
  d=Mr+beta*d;
  disp(num2str([iter rMr]));
  iter=iter+1;
end


destroy_map(tmp);

return



function[b]=make_initial_rhs(tods,tmp)
ntod=length(tods);
clear_map(tmp);
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

function[map]=apply_precon(map,wt)
mask=wt>0;
map(mask)=map(mask)./wt(mask);
return;

