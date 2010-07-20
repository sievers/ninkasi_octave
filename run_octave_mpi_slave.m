mpi_init;
while (1)
  to_eval=mpi_bcast_string();
  eval(to_eval);
end
mpi_finialize;
