function[asdf]= save_all_tod_noises(froot,tods)
myid=mpi_comm_rank+1;

if (myid==1)
  system(['mkdir ' froot ' >& /dev/null']);
end
mpi_barrier;  %make sure directory exists before proceeding.

froot(end+1)='/';
for j=1:numel(tods),
  fname=get_tod_name(tods(j));
  if isempty(fname)
    tt=num2str(j);
  else
    tt=get_tod_tags_from_names(get_tod_name(tods(j)));
  end
  save_tod_noise_banded_projvec(tods(j),[froot 'noise_' tt]);
end

