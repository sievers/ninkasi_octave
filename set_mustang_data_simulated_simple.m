function[mm]=set_mustang_data_simulated_simple(tods,mapset,input_image,varargin)
  disp(size(input_image))

sky_fac=get_keyval_default('sky_fac',1,varargin{:});
noise_fac=get_keyval_default('noise_fac',1,varargin{:});
center_map=get_keyval_default('center_map',true,varargin{:});

if ~isempty(input_image)
  if (center_map)
    input_image=center_image(input_image);
    if size(input_image,1)>min(size(mapset.skymap.map))
      nn=min(size(mapset.skymap.map));
      noff=round((size(input_image,1)-nn)/2);
      input_image=input_image(noff+1:noff+nn,noff+1:noff+nn);
    end
  end
  

  disp(size(input_image))
  disp(size(mapset.skymap.map))
  if sum(abs(size(input_image)-size(mapset.skymap.map)))==0
    disp('using original map');
    mm=input_image;
  else
    mdisp('inserting image');
    mm=insert_map_into_map(mapset.skymap.map,input_image);
  end
  mapset.skymap.map=mm;
  octave2skymap(mapset.skymap);
end

%flub=skymap2octave (mapset.skymap.mapptr);max(max(abs(flub)))


for j=1:length(tods),
  tod=tods(j);
  allocate_tod_storage(tod);
  read_tod_data(tod);
  gapfill_data_c(tod);
  array_detrend(tod);
  
  data=get_tod_data(tod);
  mat=data'*data;
  mat=0.5*(mat+mat');
  [vv,ee]=eig(mat);
  cm=data*vv(:,end);
  data_clean=data-cm*(vv(:,end)');
  vec=vv(:,end);

  assign_tod_value(tod,0);
  map2tod(mapset.skymap.mapptr,tod);

  data_sim=get_tod_data(tod);

  data_final=data_clean*noise_fac+cm*(vv(:,end)')*sky_fac+data_sim;
  free_tod_data_saved(tod);  %avoid memory leak
  set_tod_data_saved(tod,data_final);


  free_tod_storage(tod);
end

  

