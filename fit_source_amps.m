function[amps,mat]=fit_source_amps(tods,b,varargin)
if ~isfield(b,'srccat')
  error(['no srccat found in mapset in fit_source_amps.']);
end
save_tag=get_keyval_default('save_tag','map_',varargin{:});


nsrc=numel(b.srccat.amps);
mat=zeros(nsrc);
tmp_mapset=clear_mapset(b,'true');

myid=mpi_comm_rank+1;


if (myid==1)
  fid=fopen([save_tag 'srcamps_raw.dat'],'w');
  fwrite(fid,b.srccat.amps,'double');
  fclose(fid);
end


fid=fopen([save_tag '_node_' num2str(myid) '_tmp_srcmat.dat'],'w');
fid_time=fopen([save_tag '_node_' num2str(myid) '_srctimes.txt'],'w');




for j=1:nsrc,
  mdisp(['working on source ' num2str(j)]);
  tmp_mapset=clear_mapset(tmp_mapset);
  tmp_mapset.srccat.amps(j)=1;
  %we should be able to skip the MPI on the sources because we don't need to do the reduce until the end
  [tmp2,tod_times]=mapset2mapset_corrnoise_octave(tods,tmp_mapset,'skip_mpi',true,'check_empty',true,varargin{:});
  tod_times=sum(tod_times,2);
  for k=1:length(tod_times),
    fprintf(fid_time,'%10.5f ',tod_times(k));
  end
  fprintf(fid_time,'\n');
  fflush(fid_time);
  
  if isfield(tmp2,'skymap')
    destroy_map(tmp2.skymap.mapptr);
  end
  aa=tmp2.srccat.amps;
  if size(aa,1)==1,
    aa=aa';
  end
  fwrite(fid,aa,'double');
  fflush(fid);
  mat(:,j)=aa;
end
fclose(fid);
fclose(fid_time);

mat=mpi_allreduce(mat);


if isfield(tmp_mapset,'skymap')
  destroy_map(tmp_mapset.skymap.mapptr);
end


if (myid==1)
  fid=fopen([save_tag 'srcmat.dat'],'w');
  fwrite(fid,mat,'double');
  fclose(fid);
end

ii=diag(mat)~=0;

aa=b.srccat.amps;

mat_use=mat(ii,ii);
aa_use=aa(ii);
amps_small=inv(mat_use)*aa_use;
amps=0*aa;
amps(ii)=amps_small;

if (myid==1)
  fid=fopen([save_tag 'srcamps_fit.dat'],'w');
  fwrite(fid,amps,'double');
  fclose(fid);
end


