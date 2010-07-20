function[value]=mpi_reduce_map(map)

assert(class(map)=='int64');
mm=skymap2octave(map);
mm=mpi_allreduce(mm);
octave2skymap(mm,map);
   