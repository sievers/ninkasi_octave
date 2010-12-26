function[mapset,wtmap]=make_jaws_map(tods,mapset,myopts)
myid=mpi_comm_rank+1;

myopts=set_default_mapping_opts(myopts);
scale_fac=get_struct_mem(myopts,'scale_factor');
remove_common=get_struct_mem(myopts,'remove_common');
deconvolve_tau=get_struct_mem(myopts,'deconvolve_tau');
debutter=get_struct_mem(myopts,'debutter');
do_noise=get_struct_mem(myopts,'do_noise');
do_array_detrend=get_struct_mem(myopts,'detrend_array');


if (do_noise)
  noise_class=get_struct_mem(myopts,'noise_class');
  bands=get_struct_mem(myopts,'bands');
  rots=get_struct_mem(myopts,'rots');
  noise_types=get_struct_mem(myopts,'noise_types');
  assert(length(rots)==length(noise_types));
  assert(length(rots)==length(bands)-1);
  noise_scale_facs=get_struct_mem(myopts,'noise_scale_facs');
end
do_calib=get_struct_mem(myopts,'do_calib');




find_modes=get_struct_mem(myopts,'find_modes');
nbadmode=get_struct_mem(myopts,'nbadmode');
if (0) %not sure how to implement here
  find_modes_new=get_struct_mem(myopts,'find_modes_new',false);
  if (find_modes)
    nu1=get_struct_mem(myopts,'nu1');
    nu2=get_struct_mem(myopts,'nu2');
    nu3=get_struct_mem(myopts,'nu3');

    corr_root=get_struct_mem(myopts,'corr_root');
    if ~isempty(corr_root)
      ii=max(find(corr_root=='/'));
      if ~isemtpy(ii)
        corr_root=corr_root(1:ii-1);
      end
      corr_root=[corr_root '/corrs/'];
      if (myid==1)
        system(['mkdir ' corr_root]);
      end
      mpi_barrier;  %make sure directory exists before moving on
    end
    mdisp(['going to cut ' num2str(nbadmode) ' modes at ' num2str([nu1 nu2 nu3])]);
  end
end


mapset=clear_mapset(mapset);
wtmap=clear_mapset(mapset,true);
for tod_ind=1:length(tods),
  mytod=tods(tod_ind);
  allocate_tod_storage(mytod);
  read_tod_data(mytod);

  if do_calib,  calibrate_data_opts(mytod,myopts);end;
  if (scale_fac~=1), multiply_all_data(mytod,scale_fac); end


  if do_array_detrend, array_detrend(mytod);end;
  gapfill_data_c(mytod);
  if deconvolve_tau, deconvolve_tod_time_constants_c(mytod);end;
  if debutter, debutterworth_c(mytod);end;
  window_data(mytod);

  if (find_modes),
    data=get_tod_data(mytod);
    corrmat=data'*data;
    corrmat=0.5*(corrmat+corrmat');
    [vv,ee]=eig(corrmat);ee_diag=diag(ee);
    if nbadmode>0,
      nuse=nbadmode;
    else
      nuse=length(ee_diag>nbadmode^2*median(ee_diag));
    end
    

    %do noise filtering here if desired
    vv_use=vv(:,end-nuse+1:end);
    ts=data*vv_use;
    data=data-ts*vv_use';
    push_tod_data(data,mytod);
  end
  
  tod2map(mytod,mapset.skymap.mapptr);
  assign_tod_value(mytod,1);
  tod2map(mytod,wtmap.skymap.mapptr);
  free_tod_storage(mytod);
end

mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);
mapset.skymap.map=mpi_allreduce(mapset.skymap.map);



wtmap.skymap.map=skymap2octave(wtmap.skymap.mapptr);
wtmap.skymap.map=mpi_allreduce(wtmap.skymap.map);

ind=wtmap.skymap.map>0;
mapset.skymap.map(ind)=mapset.skymap.map(ind)./wtmap.skymap.map(ind);

octave2skymap(mapset.skymap);
octave2skymap(wtmap.skymap);

if (nargout==1)
  destroy_map(wtmap.skymap.mapptr);
end

