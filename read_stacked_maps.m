function[map,wt,params]=read_stacked_maps(froot,iter)

map_name=ls_fixed([froot '*_' num2str(iter) '.fits']);
if isempty(map_name)
  froot=[froot '/']
  %tt=[froot '*_' num2str(iter) '.fits']
  map_name=ls_fixed([froot '*_' num2str(iter) '.fits']);
  if isempty(map_name)
    error('unable to find any maps.');
  end
end
whos map_name
beam_name=ls_fixed([froot '*t.fits']);
wt_name=ls_fixed([froot '*s.fits']);

map_name=remove_single(map_name);
beam_name=remove_single(beam_name);
wt_name=remove_single(wt_name);

nmap=numel(map_name);

disp(['have ' num2str(numel(map_name)) ' sets, with ' num2str(numel(beam_name)) ' source maps.']);
for j=1:nmap,
  disp(map_name{j});
end

assert(numel(wt_name)==numel(map_name));
if ~isempty(beam_name)
  assert(numel(beam_name)==numel(map_name));
end

map=0;
wt=0;

for j=1:nmap,
  [params,tmp]=get_fits_projection_params(map_name{j});
  %tmp=fits_image_read(map_name{j});
  if ~isempty(beam_name)
    tmp=tmp-fits_image_read(beam_name{j});
  end
  tmp_wt=fits_image_read(wt_name{j});
  if nmap==1,
    map=tmp;
    wt=tmp_wt;
    return;
  end
  map=map+tmp.*tmp_wt;
  wt=wt+tmp_wt;
end
ii=wt>0;
map(ii)=map(ii)./wt(ii);


function[fname]=remove_single(fname)
ind=true(size(fname));
for j=1:numel(fname),
  a=strfind(fname{j},'_single_');
  if ~isempty(a)
    ind(j)=false;
  end
end
fname=fname(ind);
  