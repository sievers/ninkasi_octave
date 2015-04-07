function[mapinfo,mymap]=setup_partitioned_map(tods,bigmap,varargin)
pad=get_keyval_default('pad',10,varargin{:});

[rapix,decpix,radelt,decdelt,pv,nx,ny]=get_skymap_cea_params_c (bigmap);
mylims=get_tod_radec_lims(tods);
[mymap,mylims2]=setup_initial_map(tods,decdelt*pi/180,pad,mylims);
set_skymap_cea_simple_predef(mymap,decdelt,pv);
[rapix2,decpix2,radelt2,decdelt2,pv2,nx2,ny2]=get_skymap_cea_params_c (mymap);
i1=rapix-rapix2+1;
i2=i1+nx2-1;
j1=decpix-decpix2+1;
j2=j1-1+ny2;

mapinfo.data_patch=[i1 i2 j1 j2];
mapinfo.global_rapix=rapix;
mapinfo.global_decpix=decpix;
mapinfo.global_radelt=radelt;
mapinfo.global_decdelt=decdelt;
mapinfo.global_pv=pv;
mapinfo.global_nx=nx;
mapinfo.global_ny=ny;

nproc=mpi_comm_size;
myid=mpi_comm_rank+1;
colvec=round((0:nproc)/nproc*ny);myrows=[colvec(myid)+1 colvec(myid+1)];
mapinfo.my_patch=[1 nx myrows(1) myrows(2)];


