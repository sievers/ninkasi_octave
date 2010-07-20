addpath /home/sievers/matlab
addpath /home/sievers/act/octave

mpi_init

myrank=mpi_comm_rank;
nproc=mpi_comm_size;

targ=0;


if (1)
  mystr=num2str(myrank+1);
  str='';
  for j=0:myrank+1, 
    str=[str mystr];
  end
  nn=mpi_allreduce(1);
  str={str}
  str=repmat(str,[2*myrank 1]);
  %mystr={repmat(myrank,[1 myrank+1])};
  mystr=mpi_concatenate_cell_strings(str);
  if (myrank==0)
    disp(mystr);
   end
  %disp(['big string is .' mystr '. on ' num2str(myrank)])

end

if (0)
  x=randn(3,2+myrank);
  nn=mpi_allreduce(1);
  for j=0:nproc-1,
    if (myrank==j)
      disp(['starting ' num2str(myrank)]);
      disp(x);
    end
    pause(0.1);
  end
  y=mpi_concatenate(x,2);
  for j=0:nproc-1,
    if (myrank==j)
      disp(['ending ' num2str(myrank)]);
      disp(y);
    end
    pause(0.1);
  end
end

if (0)  
  if (myrank==targ)
    x=randn(3,4);
    disp('starting is: ')
    disp(x)
    y=mpi_bcast_array(x,targ);
    pause(0.2)
    disp('returned with:')
    disp(y);
  else
    y=mpi_bcast_array('',targ);
    pause(0.4)
    disp('received is : ')
    disp(y);
  end
end

mpi_finalize;
exit
