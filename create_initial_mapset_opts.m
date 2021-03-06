function[mapset,medians,signal_mapset,data_from_map,data_org]=create_initial_mapset_opts(tods,mapset,mapset_in,myopts)
 

myid=mpi_comm_rank+1;
nproc=mpi_comm_size;
myopts=set_default_mapping_opts(myopts);
debutter_octave=get_struct_mem(myopts,'debutter_octave');
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
dark_dirroot=get_struct_mem(myopts,'dark_dirroot','');
signal_only=get_struct_mem(myopts,'signal_only');
do_gauss=get_struct_mem(myopts,'gaussian_noise');
gauss_amp=get_struct_mem(myopts,'gauss_amp',1.0);
sim_1overf=get_struct_mem(myopts,'sim_1overf',false);
monitor_tods=get_struct_mem(myopts,'monitor_tods');
outroot=get_struct_mem(myopts,'outroot',datestr(now,30));
write_cleaned_data=get_struct_mem(myopts,'write_cleaned_data');
window_symmetric=get_struct_mem(myopts,'window_symmetric');
remove_corrnoise=get_struct_mem(myopts,'remove_corrnoise');
srccat=get_struct_mem(myopts,'srccat',[]);
do_actpol_pointing=get_struct_mem(myopts,'do_actpol_pointing',false);
free_2gamma=get_struct_mem(myopts,'free_2gamma',1);
remove_hwp=get_struct_mem(myopts,'remove_hwp',false);
read_data_new=get_struct_mem(myopts,'read_data_new',false);
gapfill_eig=get_struct_mem(myopts,'gapfill_eig',false);

abort_before_noise=get_struct_mem(myopts,'abort_before_noise',false);
skip_mpi=get_struct_mem(myopts,'skip_mpi',false);
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
  saved_noise=get_struct_mem(myopts,'saved_noise',[]);
end
do_calib=get_struct_mem(myopts,'do_calib');



add_input_map=get_struct_mem(myopts,'add_input_map');
input_scale_fac=get_struct_mem(myopts,'input_scale_fac');
restore_mapset=get_struct_mem(myopts,'restore_mapset');
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
    if (~skip_mpi)
      mpi_barrier;  %make sure directory exists before moving on
    end
  end
  mdisp(['going to cut ' num2str(nbadmode) ' modes at ' num2str([nu1 nu2 nu3])]);
end



  
mdisp(['master working on ' num2str(length(tods)) ' TODs.']);
if length(tods)==0,
   disp(['process ' num2str(myid) ' had no tods.']);
end

if do_actpol_pointing %we're going to do some precalculating of pointing.  If there's no need, then skip it.
  if ~isfield(mapset,'skymap') %if we don't couple to the sky, why do pointing?
    do_actpol_pointing=false;
  end
end

if isfield(mapset,'skymap')
  mapset.skymap=rmfield(mapset.skymap,'map');
end

for j=1:length(tods),
  tt_start=now;
  mdisp(['master working on TOD ' num2str(j) ' of ' num2str(length(tods))]);
  if length(tods)==1,
    mytod=tods;
  else
    mytod=tods(j);
  end

  if (do_actpol_pointing)

    precalc_actpol_pointing_exact(mytod,1);
    if isfield(mapset.skymap,'mapptr')
      convert_saved_pointing_to_pixellization(mytod,mapset.skymap.mapptr)
    end
    free_tod_pointing_saved(mytod,free_2gamma);
    
    if (~does_tod_have_twogamma_fit(mytod))
      precalc_actpol_pointing_exact(mytod,2);
      set_tod_twogamma_fit(mytod,'npoly_2gamma',3);
      free_tod_pointing_saved(mytod,free_2gamma);
    end

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
    if ((do_gauss)|(sim_1overf))
      if (do_gauss)        
        mdisp(['adding gaussian noise on ' get_tod_name(mytod)]);
        %crud=get_tod_data(mytod);
        %crud=randn(size(crud));
        %push_tod_data(crud,mytod);

        add_noise_to_tod_gaussian(mytod);
      end
      if (sim_1overf)
        mdisp('creating simulated 1 over f data.');
        %should now make the simulated data in-place
        make_fake_1overf_common_mode_data(mytod,myopts);
      end
    else
      mdisp('reading tod data here');
      tic;
        if (read_data_new)
          mdisp('reading data new style.');
          read_tod_data_new(mytod);
        else
          mdisp('reading data old style.');
          read_tod_data(mytod);
        end
      if (myid==1)
        toc;
      end      
      if (monitor_tods)    
        fprintf(monitor_fid,'%s\n',monitor_tag);
        fflush(monitor_fid);
      end
      
    end

    if numel(scale_fac)>1,
      fac_use=scale_fac(guess_tod_season(mytod));
      if fac_use~=1,
        multiply_all_data(mytod,fac_use);
      end      
    else
      if (scale_fac~=1)
        multiply_all_data(mytod,scale_fac);
      end
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
    

    if (gauss_amp~=1)
      dat=get_tod_data(mytod);
      dat=dat+randn(size(dat))*gauss_amp;
      push_tod_data(dat,mytod);
      clear dat;
    end


    if (hilton_noise)
      mdisp('adding hilton noise');
      array_detrend(mytod);
      gapfill_data_c(mytod);
      if debutter
        debutter_opts(mytod,myopts);
      %  if debutter_octave
      %    debutterworth_octave(mytod);
      %  else
      %    debutterworth_c(mytod);
      %  end
      end
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
      clear dat;
    end

    if (dedark)
      mdisp('dedarking');
      array_detrend(mytod);
      gapfill_data_c(mytod);
      if ~isempty(dark_dirroot)
        if iscell(dark_dirroot)
          subtract_darks_from_tod_eig(mytod,'dark_dirroot',dark_dirroot{guess_tod_season(get_tod_name(mytod))});
        else
          subtract_darks_from_tod_eig(mytod,'dark_dirroot',dark_dirroot);
        end
      else
        subtract_darks_from_tod_eig (mytod);
      end
      mdisp('finished');
    end
    
    
    if (exist('mapset_in')&(add_input_map)) | (~isempty(srccat))
      dat=get_tod_data(mytod);
      assign_tod_value(mytod,0);
      if exist('mapset_in')&(add_input_map)
        mdisp('adding input mapset into data.');
        %data_org=dat;
        mapset2tod_octave(mapset_in,mytod,j);
      end

      if ~isempty(srccat)
        inject_sources=true;
        if isfield(srccat,'inject_sources')
          if srccat.inject_sources==false,
            inject_sources=false;
          end
        end
        if inject_sources,
          mdisp('adding sources into data');
          a1=sum(sum((get_tod_data(mytod))));
          if do_actpol_pointing,
            precalc_actpol_pointing_exact(mytod);
          end
          add_srccat2tod(mytod,srccat);
          if do_actpol_pointing,          
            free_tod_pointing_saved(mytod,free_2gamma);
          end
          a2=sum(sum((get_tod_data(mytod))));
          mdisp(['finished.' num2str([a1 a1-a2])])
        end
      end
      
      if (debutter)
        if (1)
          rebutter_opts(mytod,myopts);
        else
          %mdisp('debutterworthing');
          if debutter_octave
            %mdisp('debutterworthing octave')
            debutterworth_octave(mytod,false);
          else
            %mdisp('debutterworthing c.')
            rebutterworth_c(mytod);
          end
          %mdisp('debutterworthed');
        end
      end

      if (deconvolve_tau)
        %mdisp('retauing.');
        reconvolve_tod_time_constants_c(mytod);
        %mdisp('retaued.')
      end
      %data_from_map=get_tod_data(mytod);
      %sim_dat=dat; %JLS trying to get 64 bit running
      a1=sum(sum(dat));
      dat=dat+input_scale_fac*get_tod_data(mytod);
      a2=sum(sum(dat));
      mdisp(['sums are ' num2str([a1 a2-a1 input_scale_fac])]);
      push_tod_data(dat,mytod);      
    else
      mdisp('no input mapset or source catalog');
    end
    %mdisp('doobydoobydo.');


    gapfill_data_c(mytod);
    mdisp('gapfilled.');

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
    mdisp('adding noise.');
    add_noise_to_tod(mytod);
  end

  clear dat
  clear ans
  
  dat=sum(sum(abs(get_tod_data(mytod))));mdisp(['total sum here is ' num2str(dat)]);clear dat;
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
  dat=sum(sum(abs(get_tod_data(mytod))));mdisp(['total sum here after detrend is ' num2str(dat)]);clear dat;
    
  if (debutter)
    mdisp('debutterworthing');
    if (1)
      debutter_opts(mytod,myopts);
    else
      if debutter_octave
        debutterworth_octave(mytod);
      else
        debutterworth_c(mytod);
      end
    end
  end
  if (deconvolve_tau)
    mdisp('deconvolving time constants.');
    deconvolve_tod_time_constants_c(mytod);
  end 

  if remove_hwp
    mdisp('removing half-wave plate')
    push_hwp_data=get_struct_mem(myopts,'push_hwp_data',false);
    if (push_hwp_data==false)
      warning('You have requested HWP removal, however, push_hwp_data is set to false.  You may wish to change this.');
    end
    hwp_niter=get_struct_mem(myopts,'hwp_niter',3);
    hwp_do_c=get_struct_mem(myopts,'do_hwp_az',false);
    if (hwp_do_c)
      fit_hwp_az_poly_to_data(mytod,myopts);
    else
      for iter=1:hwp_niter
        gapfill_data_c(mytod);
        [crap,crud,fitp]=fit_sines_to_hwp(mytod,myopts);
        mdisp(['on hwp iteration ' num2str(iter) ' summed fitp is ' num2str(sum(sum(abs(fitp))))]);
        clear crap
        clear crud
        clear fitp
      end
    end
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
    mdisp(['remove_corrnoise is ' num2str(remove_corrnoise)]);
    if (abort_before_noise)
      mdisp('Aborting before noise, as requested.');
      return;
    end
    
    if (gapfill_eig)
      dat=gapfill_eigenmode(get_tod_data(mytod),mytod,myopts);
      push_tod_data(dat,mytod);
    end


    use_current_data=false;
    if remove_corrnoise,
      if isfield(mapset,'corrnoise')
        mapset_tmp.corrnoise=mapset.corrnoise(j);
      end
      if isfield(mapset,'timestreams')
        mapset_tmp.timestreams=mapset.timestreams(j);
      end
      if exist('mapset_tmp')
        tic
        data_tmp=subtract_simple_guess_from_data(mytod,mapset_tmp,1);
        toc
        use_current_data=true;
      else
        warning('remove_corrnoise specified, but no valid fields found.');
        remove_corrnoise=false; %don't have anything part of the guess to pull, so why keep checking
      end
    end
    if strcmp(noise_class,'banded_projvecs')
      mdisp('setting noise to banded_projvecs');
      %if ~isempty('saved_noise')
      if numel(saved_noise)>1
        mdisp(['reading noise model from ' saved_noise]);
        read_all_tod_noises(saved_noise,mytod);
      else
        set_tod_noise_bands_projvecs(mytod,myopts);
      end
    end
    if strcmp(noise_class,'cbass_north')
      set_tod_noise_bands_cbass_north(mytod,myopts);
    end
    if strcmp(noise_class,'mustang15')
      set_tod_noise_mustang15(mytod,myopts);
    end

    if (strcmp(noise_class,'banded'))
      mdisp('setting noise to banded.');
      allocate_tod_noise_bands_c(mytod,bands);
      get_simple_banded_noise_model_c(mytod,rots,noise_types);
      %fixed, JLS, 15-mar-2011
      ndata=get_tod_ndata(mytod);
      for jj=1:length(bands)-1,
        scale_tod_band_noise_c(mytod,jj-1,ndata);
      end %end fix

      for jj=1:length(noise_scale_facs),
        scale_tod_band_noise_c(mytod,jj-1,noise_scale_facs(jj));
      end
    end      
    if strcmp(noise_class,'powlaw')
      fitp=get_mustang_noise(mytod,'order',1,'freqs',1.411,'use_current_data',use_current_data);
      set_noise_powlaw_c(mytod,fitp(1,:)/get_tod_ndata(mytod),fitp(2,:),fitp(3,:));
    end
  end
  
  if remove_corrnoise,  %gotta put back original data now.
    push_tod_data(data_tmp,mytod);
    clear data_tmp;
    clear mapset_tmp;
    clear myguess;
  end
  if ~isempty(srccat)
    mdisp('checking source restore');
    if isfield(srccat,'restore_sources')
      mdisp('flag exists');
      if srccat.restore_sources,
        mdisp('restoring sources.');
        tmpcat=srccat;
        tmpcat.amps=-1*srccat.amps;
        add_srccat2tod(mytod,tmpcat);
      end
    end
  end


  if exist('mapset_in')&(add_input_map)
    if restore_mapset
      mdisp('restoring input mapset into data.');
      dat=get_tod_data(mytod);
      assign_tod_value(mytod,0);
      mapset2tod_octave(mapset_in,mytod,j);
      if (debutter)
        if (1)
          rebutter_opts(mytod,myopts);
        else
          if debutter_octave
            debutterworth_octave(mytod,false);
          else
            rebutterworth_c(mytod);
          end
        end
      end
      if (deconvolve_tau)
        reconvolve_tod_time_constants_c(mytod);
      end
      %data_from_map=get_tod_data(mytod);
      dat=dat+get_tod_data(mytod);
    end
  end
  
  if 0
    warning('Warning!  Hand-forcing data to be signal only.')
    push_tod_data(sim_dat,mytod);
    clear sim_dat;
  end
        

  if (highpass_val>0)
    mdisp('highpassing');
    highpass_tod(mytod,highpass_val);
  end

  if (do_noise)
    mdisp('applying noise');
    %apply_banded_noise_model_c(mytod);
    apply_tod_noise_model_c(mytod);
    mdisp('applied.');
    if	check_for_nans,
      mdisp('checking for nans.');
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
  
  if window_symmetric,
    mdisp('windowing.');
    window_data(mytod);
  end

  mdisp('tod2mapset');
  mapset=tod2mapset_octave(mapset,mytod,j);
  mdisp('finished');


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
  if (do_actpol_pointing) %we cached the pointing earlier, now we need to get rid of it.
    free_tod_pointing_saved(mytod,free_2gamma); %should have already been freed, but just in case...
    free_saved_pixellization(mytod); 
  end

  if (monitor_tods)    
    fprintf(monitor_fid,'%s\n',monitor_tag);
    fflush(monitor_fid);
    fclose(monitor_fid);
  end
  
  tt_stop=now();
  mdisp(['finished processing TOD in ' num2str(86400*(tt_stop-tt_start))])
end
mdisp('master has finished his TODs');
disp(['process ' num2str(myid) ' has finished initial mapset.']);
if isfield(mapset,'skymap')
  if isfield(mapset.skymap,'partition')
    mapset.skymap=skymap2octave(mapset.skymap);
    octave2skymap(mapset.skymap);  %not sure if one needs this
  else
    %in case the map field got nuked, put it back here.
    if ~isfield(mapset.skymap,'map')
      mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);
    end
    if ~skip_mpi
      
      mapset.skymap.map=mpi_allreduce(mapset.skymap.map);
    end
    octave2skymap(mapset.skymap);
  end
end

if isfield(mapset,'srccat')
  if iscell(mapset.srccat)
    for ss=1:numel(mapset.srccat),
      if ~skip_mpi
        mapset.srccat{ss}.amps=mpi_allreduce(mapset.srccat{ss}.amps);
      end
    end
  else
    if ~skip_mpi
      mapset.srccat.amps=mpi_allreduce(mapset.srccat.amps);
    end
  end
end


if (signal_only)
  if isfield(mapset,'skyamp')
    if ~skip_mpi
      signal_mapset.skymap.map=mpi_allreduce(signal_mapset.skymap.map);
    end
    octave2skymap(signal_mapset.skymap);
  end
  if isfield(mapset,'srccat')
    if ~skip_mpi
      signal_mapset.srccat.amps=mpi_allreduce(signal_mapset.srccat.amps);
    end
  end
end

if isfield(mapset,'ground')
  if ~skip_mpi
    mpi_reduce_map(mapset.ground.groundptr);
  end
  mapset.ground.ground=skymap2octave(mapset.ground.groundptr);
end


mdisp('maps have been reduced');

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
