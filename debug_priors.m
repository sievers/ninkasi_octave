more off
addpath /cita/d/raid-sievers/sievers/matlab
mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;
format short g


%[tods_org,lims]=read_all_tod_headers(); pixsize=30/60/60*pi/180;

if (nproc==1)
  [tods_org,lims]=read_all_tod_headers(1,1); pixsize=30/60/60*pi/180;
else
  [tods_org,lims]=read_all_tod_headers('',1); pixsize=30/60/60*pi/180;
end


if (nproc==1)
  tods=tods_org(1);
  %tods=tods_org(1:4:end);
else
  tods=tods_org(myid:nproc:length(tods_org));

end

ntod=length(tods);
for j=1:ntod,
  purge_cut_detectors(tods(j));
end




map=allocate_ninkasi_skymap(pixsize,lims(1)-1e-3,lims(2)+1e-3, lims(3)-1e-3,lims(4)+1e-3);
mapset.skymap.mapptr=map;
mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);
for j=1:ntod,
  [rows,cols]=get_tod_rowcol(tods(1));
  mapset.corrnoise(j).map=randn(get_tod_ndata(tods(j)),3);
  mapset.corrnoise(j).vecs=[ones(1,how_many_dets_are_kept(tods(j)));rows'-mean(rows);cols'-mean(cols)];
end
%fid=fopen('weight_crap.out','w');fwrite(fid,size(mm,1),'int');fwrite(fid,size(mm,2),'int');fwrite(fid,mm,'double');fclose(fid);





%[mm,rhs,ww]=simple_pcg_octave_double_precon(tods(1),x);


wt=make_map_copy(map);
make_weightmap_octave(tods,wt);
weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
destroy_map(wt);
precon.skymap.map=weight;
precon.corrnoise=mapset.corrnoise;



%cmbmap=real(generate_cmb_map(size(weight),0.5));cmbmap=cmbmap';cmbmap=reshape(cmbmap,size(weight));cmbmap=150*cmbmap/mean(mean(abs(cmbmap)));
mapset2=clear_mapset(mapset,true);
%mapset2.skymap.map=cmbmap;octave2skymap(mapset2.skymap);



set_tod_noise_c(tods(1),0.1/20.0,0.00001,-1.0);
%add_noise_to_tod(tods(1));

[b,meds]=create_initial_mapset_octave(tods,mapset,'','scale_factor',1/13000);


weight_inv=weight;weight_inv(weight>0)=1./weight(weight>0);
x=clear_mapset(mapset,true);
for j=1:length(x.corrnoise),
  x.corrnoise(j).vecs=b.corrnoise(j).vecs;
  x.corrnoise(j).map(:,1)=meds{j}/mean(x.corrnoise(j).vecs(1,:));

end

x=fit_mapset_tictoc(tods,x,b,weight_inv);

return


x.skymap.map=0*x.skymap.map;
octave2skymap(x.skymap);

[b,meds]=create_initial_mapset_octave(tods,x,x,'find_modes',true,'add_noise',false,'keep_data',true);
x2=fit_mapset_tictoc(tods,x,b,weight_inv);

return










%outroot='temp_map_fft_iter_';
outroot='temp_map_debug_priors_';


%x=run_pcg_corrnoise_precon_octave(tods,x,b,precon,@apply_precon,'',20);
x=fit_mapset_tictoc(tods,x,b,weight_inv);
if (myid==1)  
  octave2skymap(x.skymap);
  write_simple_map_c(x.skymap.mapptr,[outroot 'starting_' num2str(j) '.map']);
end


  
%for j=1:0,
%  x=fit_mapset_tictoc(tods,x,b,weight_inv);
%  if (myid==1)  
%    octave2skymap(x.skymap);
%    write_simple_map_c(x.skymap.mapptr,[outroot 'starting_' num2str(j) '.map']);
%  end
%end

  
  
  
if (0)

  disp('getting priors');
  for j=1:length(tods),
    x.corrnoise(j).prior=get_smoothed_ps(x.corrnoise(j).map).^-2;
  end
  disp('got em');
  
else
  outroot=[outroot 'noprior_'];
end


for j=1:20
  [x,best]=run_pcg_corrnoise_precon_octave(tods,x,b,precon,@apply_precon,'',40);
  if (myid==1)  
    save -v7 skymap.mat x
    octave2skymap(x.skymap);
    write_simple_map_c(x.skymap.mapptr,[outroot num2str(j) '.map']);
    octave2skymap(best.skymap);
    write_simple_map_c(best.skymap.mapptr,[outroot 'best_'  num2str(j) '.map']);
    
  end
end
return

