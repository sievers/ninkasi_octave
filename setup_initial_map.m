function[map,lims]=setup_initial_map(tods,pixsize,pad)
lims=mpi_reduce_limits(get_tod_radec_lims(tods));
lims=lims+pad*pixsize*[-1 1 -1 1];
map=allocate_ninkasi_skymap(pixsize,lims(1),lims(2),lims(3),lims(4));


