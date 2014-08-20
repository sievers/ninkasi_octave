function[groundmap,precon]=setup_ground(tods,varargin)
pixsize=get_keyval_default('pixsize',30/60*pi/180,varargin{:});
pad=get_keyval_default('pad',5,varargin{:});
do_pol=get_keyval_default('do_pol',true,varargin{:});
do_precon=get_keyval_default('do_precon',true,varargin{:});

tod_map=ones(size(tods));

detlims=zeros(length(tods),4);
for j=1:length(tods),
  [dx,dy]=get_detector_offsets_actpol (tods(j));
  detlims(j,:)=[min(dx) max(dx) min(dy) max(dy)];
end
detlims=[min(detlims(:,1)) max(detlims(:,2)) min(detlims(:,3)) max(detlims(:,4))];

detlims=mpi_reduce_limits(detlims);
altazlims=get_all_tod_altaz_lims(tods);

groundlims=altazlims+detlims+pad*pixsize*[-1 1 -1 1];
groundmap=allocate_ninkasi_skymap(pixsize,groundlims(3),groundlims(4),groundlims(1),groundlims(2));
if do_pol
  set_map_polstate(groundmap,[1 1 1 0 0 0]);
end
clear_map(groundmap);


if do_precon
  precon=make_map_copy(groundmap);
  if do_pol
    set_map_polstate(precon,[1 1 1 1 1 1]);
  end
  clear_map(precon);
  for j=1:length(tods),
    allocate_tod_storage(tods(j));
    assign_tod_value(tods(j),1.0);
    tod2ground(tods(j),precon);
    free_tod_storage(tods(j));
  end
  if do_pol
    invert_pol_precon_c(precon);
  end
end


