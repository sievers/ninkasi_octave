function[mapset,vecs,cm]=create_initial_mapset_bolocam(tods,mapset,myopts,varargin)
myid=mpi_comm_rank+1;


myopts=set_default_mapping_opts(myopts);
do_noise=get_struct_mem(myopts,'do_noise');
window_symmetric=get_struct_mem(myopts,'window_symmetric');

npoly=get_struct_mem(myopts,'npoly',6);
nel_coeff=get_struct_mem(myopts,'nel_coeff',3);
nel_power=get_struct_mem(myopts,'nel_power',1);

mapset=clear_mapset(mapset);
for j=1:numel(tods),
  allocate_tod_storage(tods(j));
  read_tod_data(tods(j));
  gapfill_data_c(tods(j));
  data=get_tod_data(tods(j));

  [alt,az]=get_tod_altaz(tods(j));
  aa=alt-median(alt);aa=aa/std(aa)*std(data(:,1));
  xvec=(1:length(alt))';xvec=xvec-mean(xvec);xvec=xvec/max(abs(xvec));

  fitp=polyfit(xvec,alt,1);
  alt_resid=alt-polyval(fitp,xvec);alt_resid=alt_resid/std(alt_resid);
  alt_resid=alt_resid/std(alt_resid);


  vecs=legendre_mat(xvec,npoly);

  for kk=1:nel_power,
    avecs=alt_resid.^kk;
    for k=2:nel_coeff,
      avecs(:,end+1)=avecs(:,end).*xvec;
    end
    vecs=[vecs avecs];
  end
  rough_fitp=inv(vecs'*vecs)*(vecs'*data);

  [cm,cm_amps,vec_amps,data_clean,cm_org]=get_2pass_common_mode_wvecs(data-vecs*rough_fitp,vecs);
  assign_bad_timestreams(vecs,tods(j));
  mapset.timestreams(j).map=zeros(size(vecs,2),size(data_clean,2));
  
  if (do_noise)        
    data_tmp=data_clean;
    data_tmp(:,1)=cm;
    push_tod_data(data_tmp,tods(j));
    window_data(tods(j));
    cm_window=get_tod_data(tods(j));
    cm_window=cm_window(:,1);
    
    data_tmp=data_clean+cm*cm_amps;
    push_tod_data(data_tmp,tods(j));
    window_data(tods(j));
    set_tod_noise_bands_projvecs(tods(j),myopts);
  end
    
  push_tod_data(data_clean,tods(j));
  
  if (do_noise)
    window_data(tods(j));
    apply_tod_noise_model_c(tods(j));
    window_data(tods(j));
  end
  
  mapset=tod2mapset_octave(mapset,tods(j),j);
  free_tod_storage(tods(j));
end

if isfield(mapset,'skymap')
  mapset.skymap.map=mpi_allreduce(mapset.skymap.map);
  octave2skymap(mapset.skymap);
end






