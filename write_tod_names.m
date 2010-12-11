function[value]=write_tod_names(tods,fname)
%function[value]=write_tod_names(tods,fname)
%write the names of the tods (from get_tod_name) into fname

nproc=mpi_comm_size;
myid=mpi_comm_rank+1;

for j=length(tods):-1:1,
  tod_names(j)={get_tod_name(tods(j))};
end

if myid==1
  tag='w';
else
  tag='a';
end

for j=1:nproc,
  if myid==j
    fid=fopen(fname,tag);
    for k=1:length(tod_names),
      fprintf(fid,'%s\n',tod_names{k});
    end
    fclose(fid);
  end
  mpi_barrier;
end


