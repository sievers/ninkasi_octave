function[value]=make_reduced_weightmap_octave(tods,wt)

make_weightmap_octave(tods,wt);
weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
octave2skymap(weight,wt);
