function[weight]=make_weightmap_from_template(tods,map,varargin)
do_octave=get_keyval_default('do_octave',true,varargin{:});

wt=make_map_copy(map);
make_weightmap_octave(tods,wt);
weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
octave2skymap(weight,wt);


if (do_octave)
  %return as an octave matrix
  destroy_map(wt);
else
  %return as a map structure
  weight=wt;
end



