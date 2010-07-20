function[cm]=get_common_mode_octave(tod)
allocate_tod_storage(tod);
read_tod_data(tod);
data=get_tod_data(tod);
cm=omp_median(data,2);
free_tod_storage(tod);

