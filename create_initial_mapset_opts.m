function[mapset,medians,signal_mapset,data_from_map,data_org]=create_initial_mapset_opts(tods,mapset,mapset_in,myopts)
myid=mpi_comm_rank+1;

myopts=set_default_mapping_opts(myopts);

remove_common=get_struct_mem(myopts,'remove_common');
add_noise=get_struct_mem(myopts,'add_noise');
keep_data=get_struct_mem(myopts,'keep_data');
scale_fac=get_struct_mem(myopts,'scale_factor');
do_detrend=get_struct_mem(myopts,'detrend');
do_array_detrend=get_struct_mem(myopts,'detrend_array');
debutter=get_struct_mem(myopts,'debutter');
highpass_val=get_struct_mem(myopts,'highpass');
check_for_nans=get_struct_mem(myopts,'check_for_nans');
deconvolve_tau=get_struct_mem(myopts,'deconvolve_tau');
reverse_time=get_struct_mem(myopts,'reverse_time');
hilton_noise=get_struct_mem(myopts,'hilton_noise');
hilton_nmode=get_struct_mem(myopts,'hilton_nmode');
dedark=get_struct_mem(myopts,'dedark');
signal_only=get_struct_mem(myopts,'signal_only')
do_gauss=get_struct_mem(myopts,'gaussian_noise');
monitor_tods=get_struct_mem(myopts,'monitor_tods');
outroot=get_struct_mem(myopts,'outroot',datestr(now,30));
write_cleaned_data=get_struct_mem(myopts,'write_cleaned_data');


if (do_gauss&hilton_noise)
  warning('requested both gaussian and Hilton noise.  Choosing Gaussian.');
  hilton_noise=false;
end

if (hilton_noise)
  save_seeds=true;
else
  save_seeds=false;
end



do_noise=get_struct_mem(myopts,'do_noise');
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



add_input_map=get_struct_mem(myopts,'add_input_map');
mdisp(['add_input_map is ' num2str(add_input_map)]);



mapset=clear_mapset(mapset);
%if isfield(mapset,'corrnoise'),
%    assert(length(tods)==length(mapset.corrnoise));
%end
if exist('mapset_in')
  if isempty(mapset_in)
    signal_only=false;
    clear mapset_in;
  end
else
  signal_only=false;
end
mdisp(['signal_only is ' num2str(signal_only)])
%if isempty(mapset_in)
%  signal_only=false;
%else
  if signal_only
    signal_mapset=clear_mapset(mapset,true);
  end
%end

 

find_modes=get_struct_mem(myopts,'find_modes');
find_modes_new=get_struct_mem(myopts,'find_modes_new',false);

if (find_modes)
  nu1=get_struct_mem(myopts,'nu1');
  nu2=get_struct_mem(myopts,'nu2');
  nu3=get_struct_mem(myopts,'nu3');
  nbadmode=get_struct_mem(myopts,'nbadmode');
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



  
mdisp(['master working on ' num2str(length(tods)) ' TODs.']);
if length(tods)==0,
   disp(['process ' num2str(myid) ' had no tods.']);
end


for j=1:length(tods),
  if length(tods)==1,
    mytod=tods;
  else
    mytod=tods(j);
  end

  if (monitor_tods)
    monitor_tag=get_tod_tags_from_names(get_tod_name(mytod)) ;
    outnm=[outroot monitor_tag '.isok'];
    monitor_fid=fopen(outnm,'w');
    fprintf(monitor_fid,'%s\n',monitor_tag);
    fflush(monitor_fid);
  end

  
  allocate_tod_storage(mytod);
  if (exist('mapset_in')&(~add_input_map))
    mdisp('creating from input map.');
    mapset2tod_octave(mapset_in,mytod,j);
  else
    if (do_gauss)
      tic
      mdisp(['adding gaussian noise on ' get_tod_name(mytod)]);
      %crud=get_tod_data(mytod);
      %crud=randn(size(crud));
      %push_tod_data(crud,mytod);
      add_noise_to_tod_gaussian(mytod);
      if myid==1
        toc;
      end
    else
      mdisp('reading tod data');
      tic;
      read_tod_data(mytod);
      if (myid==1)
        toc;
      end	
    end

    if (scale_fac~=1)
      multiply_all_data(mytod,scale_fac);
    end



    if (reverse_time)
      disp('reversing time order of TOD');
      reverse_tod_uncut_regions_c(mytod);
      dat=get_tod_data(mytod);
      dat=flipud(dat);
      push_tod_data(dat,mytod);
    end
    if (do_calib)
      mdisp('calibrating');
      %calibrate_data(mytod);
      calibrate_data_opts(mytod,myopts);
    end
    
    if (hilton_noise)
      mdisp('adding hilton noise');
      array_detrend(mytod);
      gapfill_data_c(mytod);
      debutterworth_c(mytod);
      dat=get_tod_data(mytod);
      if (save_seeds)
        big_seeds(j)=round(1e9*rand);
        rand('state',big_seeds(j));
      end

      %dat=make_hilton_noise(dat,'nmode',24);
      dat=make_hilton_noise(dat,'nmode',hilton_nmode);
      if numel(dat)==1
        tt=get_tod_name(mytod);
        error(['error in hilton noise on TOD ' tt]);
      end
      push_tod_data(dat,mytod);
    end


    if (dedark)
      array_detrend(mytod);
      gapfill_data_c(mytod);      
      subtract_darks_from_tod_eig (mytod);
    end
    
    
    if exist('mapset_in')&(add_input_map)
      mdisp('adding input mapset into data.');
      dat=get_tod_data(mytod);
      %data_org=dat;
      assign_tod_value(mytod,0);
      mapset2tod_octave(mapset_in,mytod,j);
      if (debutter)
        rebutterworth_c(mytod);
      end
      if (deconvolve_tau)
        reconvolve_tod_time_constants_c(mytod);
      end
      %data_from_map=get_tod_data(mytod);
      dat=dat-get_tod_data(mytod);
      push_tod_data(dat,mytod);
      
    else
      mdisp('no input mapset');
    end


    gapfill_data_c(mytod);

    if (remove_common)
      mdisp('removing common mode');
      dat=get_tod_data(mytod);
      mns=omp_median(dat,1);
      for jjj=1:size(dat,2), dat(:,jjj)=dat(:,jjj)-mns(jjj);end;

      
      cm=omp_median(dat,2);
      for jjj=1:size(dat,2),
        dat(:,jjj)=dat(:,jjj)-cm;
      end
      push_tod_data(dat,mytod);
      clear dat;
    end


    %check for nan's after gapfilling
    if check_for_nans,
      mdisp('checking for nans');
       datamat=get_tod_data(mytod);
       nn=sum(sum(isnan(datamat)));
       clear datamat
       if nn>0
         disp(sprintf('Have %8d nans on process %3d with ndata %d on TOD %s',nn,myid,get_tod_ndata(mytod),get_tod_name(mytod)));
       end
    end
  end

  if add_noise,    
    add_noise_to_tod(mytod);
  end

  
  if (do_detrend|do_array_detrend)
    if (do_array_detrend)
      mdisp('array detrending');
      array_detrend(mytod);
      %gapfill_data_c(mytod);  %new!  
    else
      mdisp('detrending');
      detrend_data_c_c(mytod);
    end
  end
    
  if (debutter)
    mdisp('debutterworthing');
    debutterworth_c(mytod);
  end
  if (deconvolve_tau)
    mdisp('deconvolving time constants.');
    deconvolve_tod_time_constants_c(mytod);
  end
  window_data(mytod);

  
  
  
  medians{j}=tod_mean_c(mytod);
  %subtract_tod_median(mytod);
  assign_rotmat_c(mytod,0);

  if (write_cleaned_data)
    disp('writing data')
    myctime=get_tod_ctimes_from_names(get_tod_name(mytod));
    fname=sprintf('%s.%d.data',outroot,myctime);
    write_tod_data_c(mytod,fname);
  end

  if (find_modes) 
    if (find_modes_new)
      bad_modes=find_bad_modes_opts(mytod,myopts);
      [pp,r]=qr(bad_modes);
      
      mdisp('assigning rotation matrix.');
      assign_rotmat_c(mytod,pp);
      
    else
      if (1)
        corrs=get_data_correlation(mytod);
        %if (write_cleaned_data)
        %  myctime=get_tod_ctimes_from_names(get_tod_name(mytod));
        %  fname=sprintf('%s.%d.data',outroot,myctime);
        %  write_tod_data_c(mytod,fname);
        %end
        %corrs=corrs+corrs';corrs=corrs-diag(0.5*diag(corrs));
        if sum(sum(isfinite(corrs)))~= numel(corrs)
          warning(['incoming eigenvector problem on File ' get_tod_name(mytod)]);
          mapset=corrs;
          median=get_tod_data(mytod);
          return
          
        end
        
        [pp,m]=eig(corrs);m=diag(m);
        if sum(isfinite(m))~=length(m)
          error(['eigenvector problem on File ' get_tod_name(mytod)]);
        end
        mdisp('assigning rotation matrix.');
        assign_rotmat_c(mytod,pp);
      else
        a1=get_band_fft(mytod,nu1,nu2);a1=real(a1'*a1);a1=(a1+a1')/2;
        a2=get_band_fft(mytod,nu2,nu3);a2=real(a2'*a2);a2=(a2+a2')/2;
        [pp,m1,m2]=decompose_2mats(a1,a2);
        m=m1/(nu2-nu1)+m2/(nu3-nu2);  %weight per-frequency-interval,
                                      %look at total power in both matrices
      end

      for jj=1:nbadmode,
        [mm,nn]=max(m);
        badmodelist(jj)=nn;
        m(nn)=0;
      end
      bad_modes=pp(:,badmodelist);
    end
    mapset.corrnoise(j).vecs=bad_modes';
    if ~isempty(corr_root)
      nm=get_tod_name(mytod);
      while nm(end)=='/',
        nm=nm(1:end-1);
      end
      
      ii=max(find(nm=='/'));
      if ~isemtpy(ii)
        nm=nm(ii+1:end)
      end
      outname=[corr_root nm];
    end
  end

  if (do_noise)
    mdisp('doing noise');

    if (strcmp(noise_class,'banded'))
      allocate_tod_noise_bands_c(mytod,bands);
      get_simple_banded_noise_model_c(mytod,rots,noise_types);
      for jj=1:length(noise_scale_facs),
        scale_tod_band_noise_c(mytod,jj-1,noise_scale_facs(jj));
      end
    end      
    if strcmp(noise_class,'powlaw')
      fitp=get_mustang_noise(mytod,'order',1,'freqs',1.411);
      set_noise_powlaw_c(mytod,fitp(1,:)/get_tod_ndata(mytod),fitp(2,:),fitp(3,:));
    end
  end

  if (highpass_val>0)
    mdisp('highpassing');
    highpass_tod(mytod,highpass_val);
  end

  if (do_noise)
    mdisp('applying noise');
    %apply_banded_noise_model_c(mytod);
    apply_tod_noise_model_c(mytod);
    if	check_for_nans,
       datamat=get_tod_data(mytod);
       nn=sum(sum(isnan(datamat))); 

       if nn>0
         disp(['Warning - ' num2str(nn) ' extra NaNs on ' get_tod_name(mytod)]);
	 
         [rows,cols]=get_tod_rowcol(mytod);
	 disp(['bad rows are: ' num2str(rows(isnan(nn(1,:))))]);
	 disp(['bad cols are: ' num2str(rows(isnan(nn(1,:))))]);
	 return
       end
       clear datamat
     end

    %filter_tod_noise_c(mytod);
  end
  
  mdisp('tod2mapset');
  mapset=tod2mapset_octave(mapset,mytod,j);



  if (signal_only)
    if ~exist('mapset_in')
      error('signal-only requested, but no input mapset supplied.');
    end

    assign_tod_value(mytod,0);
    mapset2tod_octave(mapset_in,mytod,j);
    
    if (0)
      disp('Warning - not gapfilling.');
    else
      gapfill_data_c(mytod);
    end
    window_data(mytod);
    if (highpass_val>0)
      highpass_tod(mytod,highpass_val);
    end

    if (do_noise)
      apply_banded_noise_model_c(mytod);
    end

    if isfield(mapset,'corrnoise')
      signal_mapset.corrnoise(j).vecs=mapset.corrnoise(j).vecs;
    end

    signal_mapset=tod2mapset_octave(signal_mapset,mytod,j);
  end



  if (~keep_data)
    free_tod_storage(mytod);
  end

  if (monitor_tods)    
    fprintf(monitor_fid,'%s\n',monitor_tag);
    fflush(monitor_fid);
    fclose(monitor_fid);
  end
  
  
end

if isfield(mapset,'skymap')
  mapset.skymap.map=mpi_allreduce(mapset.skymap.map);
  octave2skymap(mapset.skymap);
end

if (signal_only)
  if isfield(mapset,'skyamp')
    signal_mapset.skymap.map=mpi_allreduce(signal_mapset.skymap.map);
    octave2skymap(signal_mapset.skymap);
  end
end


if (save_seeds)
  tod_names=cell(size(tods));
  for j=1:length(tods),
    tod_names(j)={get_tod_name(tods(j))};
  end
  tod_names=mpi_concatenate_cell_strings(tod_names);
  big_seeds=mpi_concatenate(big_seeds');
  if (myid==1)

    fid=fopen([outroot '.seeds'],'w');
    for j=1:length(tod_names),
      fprintf(fid,'%s %d %0.0f\n',tod_names{j},big_seeds(j),big_seeds(j));
    end
    fclose(fid);
  end
end




function[cm]=subtract_tod_median(tod)
dat=get_tod_data(tod);
cm=omp_median(dat,2);
dat=dat-repmat(cm,[1 size(dat,2)]);
push_tod_data(dat,tod);
return
