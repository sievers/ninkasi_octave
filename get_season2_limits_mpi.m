
more off

system('hostname');
fftw_init_threaded

addpath /home/sievers/matlab
mpi_init;
myid=mpi_comm_rank+1;

disp(['greetings from ' num2str(myid)])





nn=mpi_allreduce(1,'sum');
mdisp(['nn is ' num2str(nn)]);



do_calib=true;
highpass=0.5;



set='1';

nproc=mpi_comm_size;
format short g
%tod_names=load('season_tod_list');
%tod_names=read_text_file('equatorial_tods_trimmed.txt');


%tod_names=read_text_file('season2_ar1_tods_strip_clean2.txt');


mdisp('reading first set');
tod_names=read_text_file(['/home/sievers/act/test_priors/tods_cuts3_rising_' set  '.txt']);
mdisp('reading second set');
tod_names=[tod_names read_text_file(['/home/sievers/act/test_priors/tods_cuts3_setting_' set  '.txt'])];
mdisp(['starting with ' num2str(length(tod_names)) ' TODs.']);
tod_names=flipud(fliplr(tod_names));
for j=1:length(tod_names),
  tt=tod_names{j};if tt(end)=='/', tod_names(j)={tt(1:end-1)};end;
end

tod_names=tod_names(myid:nproc:end);
tod_names=guess_tod_name(tod_names);




decimate=0;
ntod=length(tod_names);
big_lims=zeros(ntod,4);
pixsize=30/60/60*pi/180;

tods=int64(zeros(ntod,1));




for j=1:length(tod_names),
  [tods(j),lims(j,:)]=read_all_tod_headers(tod_names(j),decimate); 
end
outroot=['season2_cuts3_set' set '_tod_limits.txt'];

if (0)
fid=fopen([outroot '_' num2str(myid)],'w');
        for j=1:ntod,
	    fprintf(fid,'%12.5f %12.5f %12.5f %12.5f %7d %s\n',lims(j,1),lims(j,2),lims(j,3),lims(j,4),get_tod_ndata(tods(j)),tod_names{j});
	end
fclose(fid);

else
for ii=1:nproc,
    if (myid==ii)
        if myid==1,
    	   fid=fopen(outroot,'w');
	else
	   fid=fopen(outroot,'a+');
        end
        for j=1:ntod,
	    fprintf(fid,'%12.5f %12.5f %12.5f %12.5f %7d %s\n',lims(j,1),lims(j,2),lims(j,3),lims(j,4),get_tod_ndata(tods(j)),tod_names{j});
	end
    fclose(fid);
    end
    assert(mpi_allreduce(1)==nproc);
end
endif
mpi_finalize;
return

