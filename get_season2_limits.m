addpath /home/sievers/matlab/
mpi_init
more off
tod_names=read_text_file('season_tod_list');
%tod_names=tod_names(1:100);
tod_names=guess_tod_name(tod_names);
decimate=6;


fid=fopen('season2_ar1_limits.txt','w');
n=length(tod_names)
big_lims=zeros(n,4);
for j=1:length(tod_names)
  [tod,big_lims(j,:)]=read_tod_header(tod_names{j},'',decimate);
  fprintf(fid,'%12.4f %12.4f %12.4f %12.4f %s\n',big_lims(j,1),big_lims(j,2),big_lims(j,3),big_lims(j,4),tod_names{j});
  fflush(fid);
end
fclose(fid);

