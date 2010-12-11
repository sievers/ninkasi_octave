function[value]=write_char_cells_mpi(cells,fname)
big_cells=mpi_concatenate_cell_strings(cells);
myid=mpi_comm_rank;
if myid==0,
  fid=fopen(fname,'w');
  for j=1:length(big_cells),
    fprintf(fid,'%s\n',big_cells{j});
  end
  fclose(fid);
end
