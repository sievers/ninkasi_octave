aaa=whos;
existing_vars={};
for j=1:length(aaa),
  myvar=aaa(j);
  if (strcmp(myvar.name,'ans')==0)&(strcmp(myvar.name,'aaa')==0)
    if eval(['ischar(' myvar.name ')'])
      nextline=[myvar.name ' = ''' eval(myvar.name) ''''];
    else
      nextline=[myvar.name ' = ' eval(['num2str(' myvar.name ');'])];
    end
    existing_vars(end+1)=nextline;
  end
end

                


#addpath /home/sievers/matlab

more off
crash_dumps_octave_core(false);



fftw_init_threaded

mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;

if nproc>1
  ignore_function_time_stamp ("all")
end
