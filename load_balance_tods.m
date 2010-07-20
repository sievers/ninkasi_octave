function[tod_names]=load_balance_tods(tod_names,opts)
ncpu=mpi_comm_size;
myid=mpi_comm_rank+1;

tod_names=guess_tod_name(tod_names);
[tods,lims]=read_all_tod_headers(tod_names,opts.decimate); 
tods=setup_tods_opts(tods,opts);

if isempty(tods)
  mywork=[];
else
  mywork(length(tods),1)=0;
end

for j=1:length(tods),
  mywork(j,1)=get_tod_ndata(tods(j))*length(get_tod_rowcol(tods(j)));
  assert(mywork(j)>0);  
end

if isempty(tods)
  tod_names={};
else
  tod_names=get_tod_names(tods);
end

big_work=mpi_concatenate(mywork);
big_names=mpi_concatenate_cell_strings(tod_names);
total_work=sum(big_work);

cpu_work=mpi_concatenate(sum(mywork));
if myid==1
  mdisp(['pre balancing work min mean max (gigapixels): ' num2str([min(cpu_work) mean(cpu_work) max(cpu_work)]/1e9)]);
end




assert(total_work==mpi_allreduce(sum(mywork)));  %make sure that ended up OK
%disp(['total work is ' num2str(total_work) ' vs ' num2str(mpi_allreduce(sum(mywork))) ' with ' num2str([length(big_names) length(big_work)])]);




assert(length(big_names)==length(big_work));
[big_work,ind]=sortrows(big_work);

big_names=big_names(ind);

if size(big_names,1)==1,
  big_names=big_names';
end
nn=floor(length(big_work)/ncpu);
big_names=flipud(big_names);
big_work=flipud(big_work);

for j=1:2:nn,
  ii=(j-1)*ncpu+1:(j*ncpu);
  big_names(ii)=flipud(big_names(ii));
  big_work(ii)=flipud(big_work(ii));
end

tod_names=big_names(myid:ncpu:end);
mywork=big_work(myid:ncpu:end);
cpu_work=mpi_concatenate(sum(mywork));
if myid==1
  mdisp(['work min mean max (gigapixels): ' num2str([min(cpu_work) mean(cpu_work) max(cpu_work)]/1e9)]);
end
assert(sum(cpu_work)==total_work);

disp([num2str(myid) ' is working on ' num2str(length(tod_names)) ' tods']);
%crud=mpi_concatenate_cell_strings(tod_names);
if myid==0,
  fid=fopen('test_names.txt','w');
  for j=1:length(crud),
    fprintf(fid,'%s\n',crud{j});
  end
  fclose(fid);
end


