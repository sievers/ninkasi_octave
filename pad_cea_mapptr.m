function[value]=pad_cea_mapptr(mapptr)
fwee=size(skymap2octave(mapptr));
myvec=find_all_fft_pads_v2(max(fwee));
new_nx=min(myvec(myvec>=fwee(1)));
new_ny=min(myvec(myvec>=fwee(2)));
resize_skymap_cea_c(mapptr,new_nx,new_ny);

