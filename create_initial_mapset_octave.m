function[mapset,medians,dat_org,dat_from_map]=create_initial_mapset_octave(tods,mapset,mapset_in,varargin)
myid=mpi_comm_rank+1;



add_noise=get_keyval_default('add_noise',false,varargin{:});  %please assign the noise levels before coming in!
keep_data=get_keyval_default('keep_data',false,varargin{:});
scale_fac=get_keyval_default('scale_factor',1,varargin{:});
do_detrend=get_keyval_default('detrend',false,varargin{:});
do_array_detrend=get_keyval_default('detrend_array',false,varargin{:});
debutter=get_keyval_default('debutter',false,varargin{:});
highpass_val=get_keyval_default('highpass',0,varargin{:});
check_for_nans=get_keyval_default('check_for_nans',true,varargin{:});
deconvolve_tau=get_keyval_default('deconvolve_tau',false,varargin{:});
reverse_time=get_keyval_default('reverse_time',false,varargin{:});
hilton_noise=get_keyval_default('hilton_noise',false,varargin{:});
dedark=get_keyval_default('dedark',false,varargin{:});


do_noise=get_keyval_default('do_noise',false,varargin{:});
if (do_noise)
  bands=get_keyval_default('bands',[-1 3 300],varargin{:});
  rots=get_keyval_default('rots',[0 0],varargin{:});
  noise_types=get_keyval_default('noise_types',[1 0],varargin{:});
  assert(length(rots)==length(noise_types));
  assert(length(rots)==length(bands)-1);
  noise_scale_facs=get_keyval_default('noise_scale_facs',[],varargin{:});
end
do_calib=get_keyval_default('do_calib',false,varargin{:});



add_input_map=get_keyval_default('add_input_map',false,varargin{:});
mdisp(['add_input_map is ' num2str(add_input_map)]);



mapset=clear_mapset(mapset);
%if isfield(mapset,'corrnoise'),
%    assert(length(tods)==length(mapset.corrnoise));
%end
if exist('mapset_in')
  if isempty(mapset_in)
    clear mapset_in;
  end
end

find_modes=get_keyval_default('find_modes',false,varargin{:});
if (find_modes)
  nu1=get_keyval_default('nu1',0.01,varargin{:});
  nu2=get_keyval_default('nu1',4,varargin{:});
  nu3=get_keyval_default('nu1',30,varargin{:});
  nbadmode=get_keyval_default('nbadmode',3,varargin{:});
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
  
  allocate_tod_storage(mytod);
  if (exist('mapset_in')&(~add_input_map))
    mdisp('creating from input map.');
    mapset2tod_octave(mapset_in,mytod,j);
  else
    mdisp('reading tod data');
    tic;
    read_tod_data(mytod);
    if (myid==1)
      toc;
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
      calibrate_data(mytod);
    end
    
    if (hilton_noise)
      array_detrend(mytod);
      gapfill_data_c(mytod);
      debutterworth_c(mytod);
      dat=get_tod_data(mytod);
      dat=make_hilton_noise(dat,'nmode',24);
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
      dat_org=dat;
      assign_tod_value(mytod,0);
      mapset2tod_octave(mapset_in,mytod,j);
      if (debutter)
        rebutterworth_c(mytod);
      end
      if (deconvolve_tau)
        reconvolve_tod_time_constants_c(mytod);
      end
      dat_from_map=get_tod_data(mytod);
      dat=dat+get_tod_data(mytod);
      push_tod_data(dat,mytod);
    else
      mdisp('no input mapset');
    end


    gapfill_data_c(mytod);

    %check for nan's after gapfilling
    if check_for_nans,
       datamat=get_tod_data(mytod);
       nn=sum(sum(isnan(datamat)));
       clear datamat
       disp(sprintf('Have %8d nans on process %3d with ndata %d on TOD %s',nn,myid,get_tod_ndata(mytod),get_tod_name(mytod)));
    end
  end

  if add_noise,    
    add_noise_to_tod(mytod);
  end

  
  if (do_detrend|do_array_detrend)
    if (do_array_detrend)
      mdisp('array detrending');
      array_detrend(mytod);
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
  if (find_modes)    
    if (1)
      corrs=get_data_correlation(mytod);
      %corrs=corrs+corrs';corrs=corrs-diag(0.5*diag(corrs));
      [pp,m]=eig(corrs);m=diag(m);
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
    mapset.corrnoise(j).vecs=bad_modes';
  end

  if (do_noise)
    %allocate_tod_noise_bands_c(mytod,[-1  3 300]);
    %get_simple_banded_noise_model_c(mytod,[ 1 0],[  1 0]);
    allocate_tod_noise_bands_c(mytod,bands);
    get_simple_banded_noise_model_c(mytod,rots,noise_types);
    for jj=1:length(noise_scale_facs),
      scale_tod_band_noise_c(mytod,jj-1,noise_scale_facs(jj));
    end
							 
    %fit_tod_noise_c(mytod);
  end

  
  if (highpass_val>0)
    highpass_tod(mytod,highpass_val);
  end

  if (do_noise)
    apply_banded_noise_model_c(mytod);
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
  
  mapset=tod2mapset_octave(mapset,mytod,j);
  if (~keep_data)
    free_tod_storage(mytod);
  end
  
end

mapset.skymap.map=mpi_allreduce(mapset.skymap.map);
octave2skymap(mapset.skymap);



function[cm]=subtract_tod_median(tod)
dat=get_tod_data(tod);
cm=omp_median(dat,2);
dat=dat-repmat(cm,[1 size(dat,2)]);
push_tod_data(dat,tod);
return
