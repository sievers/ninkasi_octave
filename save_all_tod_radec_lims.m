function[value]=save_all_tod_radec_lims(tods,outroot)
%save the tod radec limits to a file.
[fwee,lims]=get_tod_radec_lims(tods);

myid=mpi_comm_rank+1;
nproc=mpi_comm_size;
outname=[outroot 'radecs.txt'];
for j=1:nproc,
  if myid==j,
    if myid==1
      fid=fopen(outname,'w');
    else
      fid=fopen(outname,'a');
    end
    for k=1:numel(tods),
      fprintf(fid,'%s  %14.6f %14.6f %14.6f %14.6f\n',get_tod_name(tods(k)),lims(k,1),lims(k,2),lims(k,3),lims(k,4));
    end
    fclose(fid);
  end
  mpi_barrier;
end

