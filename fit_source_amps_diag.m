function[amps,curve]=fit_source_amps_diag(tods,b,varargin)
if ~isfield(b,'srccat')
  error(['no srccat found in mapset in fit_source_amps.']);
end
if iscell(b.srccat),
  save_tag=get_keyval_default('save_tag','map_',varargin{:});
  amps=cell(size(b.srccat));
  curve=cell(size(b.srccat));
  for j=1:numel(b.srccat),
    mytag=[save_tag '_set_' num2str(j)];
    bb=b;
    bb.srccat=b.srccat{j};
    [fwee,fum]=fit_source_amps_diag(tods,bb,'save_tag',mytag,varargin{:});
    amps(j)={fwee};
    curve(j)={fum};
  end
  return
end



save_tag=get_keyval_default('save_tag','map_',varargin{:});


nsrc=numel(b.srccat.amps);

tmp_mapset=clear_mapset(b,'true');

myid=mpi_comm_rank+1;


if (myid==1)
  fid=fopen([save_tag 'srcamps_raw.dat'],'w');
  fwrite(fid,b.srccat.amps,'double');
  fclose(fid);
end



curve=zeros(nsrc,1);

srccat=b.srccat;
  
aa=now;
for ii=1:length(tods),
  mytod=tods(ii);
  do_i_hit=tod_hits_srccat(mytod,srccat);
  mdisp(['working on tod' num2str(ii) ' which hits ' num2str(sum(do_i_hit)) ' sources.']);



  for j=1:nsrc,
    if do_i_hit(j),
      %mdisp(['working on source ' num2str(j)])
      tmp_srccat=srccat;
      tmp_srccat.ra=srccat.ra(j);
      tmp_srccat.dec=srccat.dec(j);
      tmp_srccat.amps=1;
      tmp_mapset.srccat=tmp_srccat;
      [tmp2,tod_times]=mapset2mapset_corrnoise_octave(mytod,tmp_mapset,'skip_mpi',true,'check_empty',false,varargin{:});


      curve(j)=curve(j)+tmp2.srccat.amps;
    end
  end
end
bb=now;
disp(['node ' num2str(myid) ' took ' num2str(86400*(bb-aa)) ' seconds to do sources.']);
curve=mpi_allreduce(curve);

if size(b.srccat.amps,1)==size(curve,1)
  amps=b.srccat.amps./curve;
else
  amps=b.srccat.amps'./(curve);
end
amps(curve==0)=0;


if (myid==1)
  fid=fopen([save_tag '_srccurve.dat'],'w');
  fwrite(fid,curve,'double');
  fclose(fid);
  
  fid=fopen([save_tag '_srcamps_fit_diag.dat'],'w');
  for j=1:nsrc,
    err=0;
    if curve(j)>0,
      err=1/sqrt(curve(j));
    end
    fprintf(fid,'%14.5g %14.5g\n',amps(j),err);
  end
  fclose(fid);
end

