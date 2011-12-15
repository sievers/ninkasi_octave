function[asdf]=write_tod_lims(tod_names,pfile,outfile)

ntod=numel(tod_names);


[tods,lims,all_lims]=read_all_tod_headers_vararg(tod_names,'pointing',pfile);
ii=all_lims(:,2)<all_lims(:,1);
all_lims(ii,2)=all_lims(ii,2)+2*pi;
mystr=cell(ntod,1);
for j=1:ntod,
  ll=get_tod_name(tods(j));
  ll=[ll '  ' sprintf(' %14.6f ',all_lims(j,:))];
  mystr(j)={ll};
end

mpi_barrier;
mdisp('writing cells')
write_char_cells_mpi(mystr,outfile);
mpi_barrier;
mdisp('wrote cells');