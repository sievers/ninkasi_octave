function[weight]=simple_make_weightmap_octave(tods,mapset,outroot)

wt=make_map_copy(mapset.skymap.mapptr);
make_weightmap_octave(tods,wt);



weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
mdisp('made weightmap');
octave2skymap(weight,wt);
if (mpi_comm_rank==0)
  write_map(wt,[outroot 'weights']);
end

destroy_map(wt);

