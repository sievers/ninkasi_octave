function[value]=write_tod_dists(tods,ra,dec,fname)

for j=length(tods):-1:1,
  dists(j)={get_tod_dist_vec_c(ra,dec,tods(j))};
end

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
      fprintf(fid,'%s ',tod_names{k});
      dd=dists{k};
      for kk=1:length(dd),
        fprintf(fid,'%12.7f ',dd(kk));
      end
      fprintf(fid,'%12.7f\n',sum(dd));
    end
    fclose(fid);
  end
  mpi_barrier;
end

