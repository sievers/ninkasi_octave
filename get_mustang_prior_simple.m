function[fitp,datft,x,noise]=get_mustang_prior_simple(tod,corrvecs,varargin)
allocate_tod_storage(tod);
read_tod_data (tod);
gapfill_data_c (tod);
array_detrend (tod);
data=get_tod_data(tod);


mapset_in=get_keyval_default('mapset',[],varargin{:});
if ~isempty(mapset_in)
  mapset_use.skymap=mapset_in.skymap;
  assign_tod_value(tod,0.0);
  octave2skymap(mapset_use.skymap);  
  mapset2tod_octave(mapset_use,tod,1);
  data=data-get_tod_data(tod);
  clear mapset_use;
  clear mapset_in;
end


free_tod_storage(tod);


ts=pull_bad_timestreams(tod);


dd=remove_timestreams_from_corrnoise(data,ts);
my_vecs=dd*corrvecs';
global x;


[datft,x]=get_corrnoise_ft_nu(my_vecs);
mypow=1.0;
datft=datft.^mypow;


slope=-1.8;
vecs=[ones(size(x)) x.^slope];
for j=size(datft,2):-1:1,
  pp(:,j)=linfit(datft(:,j),vecs);
end
for iter=1:5,
  for j=size(datft,2):-1:1,
    noise=(pp(1,j)+pp(2,j)*x.^slope);
    pp(:,j)=linfit(datft(20:end,j),vecs(20:end,:));
  end
end



fitp=[pp(1,:); ones(1,size(datft,2))*slope;pp(2,:)];

global rescale





global datasqr;
for j=1:size(datft,2),
  datasqr=datft(:,j);

  rescale=pp(:,j);
  fitp(:,j)=fminunc(@get_like_falpha,fitp(:,j));  
  fitp(1,j)=fitp(1,j)*rescale(1);
  fitp(3,j)=fitp(3,j)*rescale(2);
  rescale=[fitp(1,j);fitp(3,j)];
  fitp(1,j)=1;
  fitp(3,j)=1;
  fitp(:,j)=fminunc(@get_like_falpha,fitp(:,j));  
  fitp(1,j)=fitp(1,j)*rescale(1);
  fitp(3,j)=fitp(3,j)*rescale(2);

end


return


function[like]=get_like_falpha(fitp)
global rescale
global x
global datasqr

model=fitp(1)*rescale(1)+fitp(3)*rescale(2)*(x.^fitp(2));

like=sum(datasqr./model)+sum(log(model));

if ~isreal(like)
  like=1e10+abs(like);
end

%disp(sprintf('%14.2f %12.4e %12.4g %12.4f',like,amp,knee,slope));


