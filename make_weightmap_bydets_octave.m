function[value]=make_weightmap_bydets_octave(tods,map,do_window,varargin)

do_new_pointing=get_keyval_default('do_new_pointing',false,varargin{:})
do_reduce=get_keyval_default('mpi_reduce',false,varargin{:});


myrows=get_keyval_default('rows',[],varargin{:});
mycols=get_keyval_default('cols',[],varargin{:});
myweights=get_keyval_default('det_weights',[],varargin{:});

clear_map(map);

if isempty('do_window')
  clear do_window;
end


if ~exist('do_window')
  do_window=true;
end

for j=1:length(tods),    
  mytod=tods(j);
  allocate_tod_storage(mytod);
  [rr,cc]=get_tod_rowcol(mytod);
  ww=0*rr;
  for jj=1:length(rr),
    myind=find( (rr(jj)==myrows)&(cc(jj)==mycols));
    assert(numel(myind)==1);  %if this fails, you didn't pass in a well-behaved set of row/columns
    ww(jj)=myweights(myind);
  end
  assign_tod_detector_values_c(tods,ww);;
  
  if do_window,
    window_data(mytod);
    window_data(mytod);
  end
  
  
  if (do_new_pointing)
    mdisp('doing pointing')
    precalc_actpol_pointing_exact(mytod);
    convert_saved_pointing_to_pixellization(mytod,map);
    free_tod_pointing_saved(mytod);
  end
  
  mdisp('assigned value');
  if (is_map_polarized(map))
    tod2polmap(mytod,map);
  else
    tod2map(mytod,map);
  end
  mdisp('projected');
  
  if (do_new_pointing)
    %free_tod_pointing_saved(mytod);
    free_saved_pixellization(mytod);
    
  end


  free_tod_storage(mytod);    
  %disp('freed');
end

if (do_reduce)
  mdisp('reducing')
  weight=skymap2octave(map);
  weight=mpi_allreduce(weight);
  octave2skymap(weight,map);
  clear weight;
end

  
