function[map,old_mapset,pnb,map2,bnp]=fit_mapset_tictoc(tods,map,b,weight_inv)

%make a copy of the incoming map, with the correlated noise zeroed out.
map2=clear_mapset(map,'true');
map2.skymap.map=map.skymap.map;
octave2skymap(map2.skymap);

bnp=mapset2mapset_corrnoise_octave(tods,map2);



%first, solve for the current correlated noise estimates.
new_mapset=clear_mapset(map,'true');
for j=1:length(tods),
  %little_mat=map.corrnoise(j).vecs*map.corrnoise(j).vecs';
  %new_mapset.corrnoise(j).map=(b.corrnoise(j).map- bnp.corrnoise(j).map)*inv(little_mat);
  if isfield(map.corrnoise(j),'prior')
    new_mapset.corrnoise(j).map=fit_corrnoise_wpriors(map.corrnoise(j).vecs,(b.corrnoise(j).map- bnp.corrnoise(j).map),map.corrnoise(j).prior);
  else
    new_mapset.corrnoise(j).map=fit_corrnoise_wpriors(map.corrnoise(j).vecs,(b.corrnoise(j).map- bnp.corrnoise(j).map));
  end
  
end

  

%now solve for the current skymap, given current correlations.  new_mapset should have sky=0 at this point
pnb=mapset2mapset_corrnoise_octave(tods,new_mapset);


old_mapset=new_mapset;


rhs=b.skymap.map-pnb.skymap.map;
new_mapset.skymap.map=rhs.*weight_inv;


%if (mapsetdotmapset(map2,bnp)>0)
if (0)
  Arhs=mapset2mapset_corrnoise_octave(tods,new_mapset);
  % 1)=map2/bnp.  2)=old_mapset/pnb  3)=new_mapset/Arhs
  mat(1,1)=mapsetdotmapset(map2,bnp);
  mat(1,2)=mapsetdotmapset(map2,pnb);
  mat(1,3)=mapsetdotmapset(map2,Arhs);
  mat(2,2)=mapsetdotmapset(old_mapset,pnb);
  mat(2,3)=mapsetdotmapset(old_mapset,Arhs);
  mat(3,3)=mapsetdotmapset(new_mapset,Arhs);
  mat(2,1)=mat(1,2);
  mat(3,1)=mat(1,3);
  mat(3,2)=mat(2,3);
  assert(mat==mat');
  vec(1,1)=mapsetdotmapset(map2,b);
  vec(2,1)=mapsetdotmapset(old_mapset,b);
  vec(3,1)=mapsetdotmapset(new_mapset,b);
  disp(mat);
  disp(vec);
  val=inv(mat)*vec;  
  new_mapset=mapset_axpy(new_mapset,old_mapset,val(2),val(3));
  new_mapset=mapset_axpy(new_mapset,map2,val(1));
end






destroy_map(map2.skymap.mapptr);
map.skymap.map=new_mapset.skymap.map;
map.corrnoise=new_mapset.corrnoise;
octave2skymap(map.skymap);
destroy_map(new_mapset.skymap.mapptr);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555


function[fit]=fit_corrnoise_vecs_fft(vecs,map,prior)
little_mat=vecs*vecs';
mapft=map;
for j=1:size(map,2),
  mapft(:,j)=fft(mapft(:,j));
end
if (~exist('prior'))
  fitft=mapft*inv(little_mat);
else
  
  val1=multiply_corrnoise_wprior_c(real(mapft)',1./(prior.^2)',little_mat);
  val2=multiply_corrnoise_wprior_c(imag(mapft)',1./(prior.^2)',little_mat);
  fitft=val1'+i*val2';
  %for j=1:length(map),
  %  fitft(j,:)=mapft(j,:)*inv(little_mat+diag(1./prior(j,:).^2));
  %end
end

fit=zeros(size(fitft));
for j=1:size(fitft,2),
  fit(:,j)=real(ifft(fitft(:,j)));
end

  

