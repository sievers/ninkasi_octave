function[x_accum]=run_pcg_multiscale(tods,x,b,precon,apply_precon,apply_prior_to_mapset,varargin)
myid=mpi_comm_rank+1;

if numel(varargin)>1,
  clear myopts;
  for j=1:2:numel(varargin)
    eval(['myopts.' varargin{j} ' = varargin{j+1};']);
  end
else
  myopts=varargin{1};
end

if isstruct(precon)
  if isfield(precon,'skymap')
    if ~isfield(precon.skymap,'mapptr')
      if isfield(b,'skymap')
        precon.skymap.mapptr=b.skymap.mapptr;
      end
    end
  end
end

scales=get_struct_mem(myopts,'scales',[0 1 2 3 4 0]);
iters_default=50;
iters=get_struct_mem(myopts,'iters',iters_default);
if numel(iters)<numel(scales),
  iters(end+1:numel(scales))=iters(end);
end
disp([scales;iters]);

global weight_template  %let's not worry about scaling the weight template for now
weight_template=[];


%b_cur=b;
%x_cur=x;
x_accum=clear_mapset(x);

for curscale=1:numel(scales),
  pcg_opts=myopts;
  pcg_opts.save_tag=[myopts.save_tag 'scale_' num2str(curscale) '_'];
  pcg_opts.maxiter=iters(curscale);
  myres=scales(curscale);
  b_use=downres_mapset(b,myres);
  x_use=downres_mapset(x,myres,true); %true is to make sure we take the mean, not the sum
  precon_use=downres_mapset(precon,myres);
  
  [x_use]=run_pcg_corrnoise_precon_octave(tods,x_use,b_use,precon_use,apply_precon,apply_prior_to_mapset,pcg_opts);
  if (myid==1)
    if isfield(x_use,'skymap')
      write_map(x_use.skymap.mapptr,[pcg_opts.save_tag 'ending']);
    end
  end

  pcg_opts.new_mapptr=b.skymap.mapptr;
  tmp_mapset=mapset2mapset_corrnoise_octave(tods,x_use,pcg_opts);
  b=add_mapset(b,tmp_mapset,-1);
  if (isfield(tmp_mapset,'skymap'))
    destroy_map(tmp_mapset.skymap.mapptr);
  end
  clear tmp_mapset



  if isfield(x_use,'skymap')
    disp('upresing')
    tmp=x_use.skymap.mapptr;
    octave2skymap(x_use.skymap);
    x_use.skymap.mapptr=upres_mapptr(x_use.skymap.mapptr,myres);
    x_use.skymap.map=skymap2octave(x_use.skymap.mapptr);
    destroy_map(tmp);
  end
  add_mapset(x_accum,x_use,1);
  if (isfield(x_accum,'skymap'))
    if (myid==1)
      octave2skymap(x_accum.skymap);
      write_map(x_accum.skymap.mapptr,[pcg_opts.save_tag 'accum']);
    end
    destroy_map(x_use.skymap.mapptr);
    destroy_map(b_use.skymap.mapptr);
    
  end  
  x=clear_mapset(x);
  clear x_use;
  clear b_use;


end



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[mapset]=downres_mapset(mapset,res,rescale)
if ~isstruct(mapset)
  return
end

if ~isfield(mapset,'skymap')
  return
end
if ~exist('rescale')
  rescale=false;
end


mapset.skymap.mapptr=make_map_copy(mapset.skymap.mapptr);
octave2skymap(mapset.skymap);
scale_fac=1;
for j=1:res,
  map_deres=deres_map(mapset.skymap.mapptr);
  if (rescale)
    scale_fac=scale_fac/0.25;
  end
  destroy_map(mapset.skymap.mapptr);
  mapset.skymap.mapptr=map_deres;
end
mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);
if (scale_fac ~=1)
  mapset.skymap.map=mapset.skymap.map*scale_fac;
  octave2skymap(mapset.skymap);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[mapptr]=upres_mapptr(mapptr,res)
mapptr=make_map_copy(mapptr);
for j=1:res,
  tmp=upres_map(mapptr);
  destroy_map(mapptr);
  mapptr=tmp;
end

