function[x]=run_pcg_corrnoise_precon_octave(tods,x,b,precon,fun,priorfun,varargin)
%function[x,best]=run_pcg_corrnoise_precon_octave(tods,x,b,precon,fun,priorfun,varargin)

%x is a mapset, could be clear
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;

maxiter=get_keyval_default('maxiter',50,varargin{:});
tol=get_keyval_default('tol',1e-8,varargin{:});
save_interval=get_keyval_default('save_interval',1e-8,varargin{:});
save_tag=get_keyval_default('save_tag','map_',varargin{:});
dad_thresh=get_keyval_default('dad_thresh',1e3,varargin{:});
write_times=get_keyval_default('write_times',false,varargin{:});
profile=get_keyval_default('profile',false,varargin{:});

if (myid==1)
  disp(['length of varargin is ' num2str(length(varargin))]);
  disp(['save interval is ' num2str(save_interval)]);
  disp(['save tag is ' save_tag]);
  disp(['tol is ' num2str(tol)]);
  disp(['maxiter is ' num2str(maxiter)]);
end


if exist('priorfun')
  if isempty(priorfun)
    clear priorfun;
  end
end
if exist('precon')
  if isempty(precon)
    clear precon;
  end
end



if exist('precon')
  do_precon=true;
else
  do_precon=false;
end

if ~exist('b')
  b='';
end
if isempty(b)
  b=create_initial_mapset_octave(tods,x);
end



[ax,tod_times,node_times]=mapset2mapset_corrnoise_octave(tods,x,varargin{:});
if myid==1
  if isfield(ax,'skymap')
    write_map(ax.skymap,[save_tag '_ax']);
  end
  fieldnames(ax)
end


if (write_times)
  for j=1:mpi_comm_size,
    if (myid==j)
      if (myid==1)
        fid=fopen([save_tag '.times'],'w');
        fid_node=fopen([save_tag '.node_times'],'w');

      else
        fid=fopen([save_tag '.times'],'a');
        fid_node=fopen([save_tag '.node_times'],'a');
      end
      for k=1:length(tods),
        tod_name=get_tod_name(tods(k));
        [rows,cols]=get_tod_rowcol(tods(k));
        %fprintf(fid,'%s %14.4f %6d %8d %s\n',tod_name,tod_times(k),length(rows),get_tod_ndata(tods(k)),getenv('HOSTNAME'));
        fprintf(fid,'%s %14.4f %14.4f %14.4f %14.4f %6d %8d %4d\n',tod_name,tod_times(k,1),tod_times(k,2),tod_times(k,3),tod_times(k,4),length(rows),get_tod_ndata(tods(k)),myid);
      end
      fclose(fid);
      fprintf(fid_node,'%14.4f %4d\n',node_times,myid);
      fclose(fid_node);
      %system('env');
    end
    mpi_barrier;
  end
end


if exist('priorfun')
  %disp('evaluating prior');  
  ax=feval(priorfun,ax,x);
  %disp('evaluated.');
end

r=add_mapset(b,ax,-1);
if do_precon,
  %disp('evaluating precon');
  Mr=feval(fun,r,precon);
  %disp('evaluated');
else
  Mr=r;
end

d=Mr;



rMr=mapsetdotmapset(r,d);
r0sqr=rMr;
iter=1;

if (myid==1)
  disp(['r0sqr is ' num2str(r0sqr)]);
end

best_rMr=r0sqr;



last_failed=0;


if (profile)
  logid=fopen([save_tag '.profile.node_' num2str(myid)],'w');
  for j=1:length(tods),
    fprintf(logid,'%s\n',get_tod_name(tods(j)));
  end
  fflush(logid);
end


while ((rMr>r0sqr*tol)&(iter<maxiter)),
  %tic;
  if (profile)  %if profiling, make sure we all start at the same time
    mpi_barrier;  
  end
    aa=now;
    %sz1=size(d.skymap.map);
    [Ad,tod_times,node_time,mpi_time]=mapset2mapset_corrnoise_octave(tods,d,varargin{:});
    cc=now;

    slowest_time=mpi_allreduce(node_time,'max');
    if node_time==slowest_time
      bad_id=myid;
    else
      bad_id=0;
    end
    bad_node=mpi_allreduce(bad_id);
    total_time=mpi_allreduce(node_time);

    %sz2=size(Ad.skymap.map); 
    %if min(sz1==sz2)==0, error(['spot 1 size mismatch: ' num2str(sz1) ' vs ' num2str(sz2)]);end;
    if exist('priorfun')
      Ad=feval(priorfun,Ad,d);
    end
    %sz2=size(Ad.skymap.map); 
    %if min(sz1==sz2)==0, error(['spot 2 size mismatch: ' num2str(sz1) ' vs ' num2str(sz2)]);end;

    dAd=mapsetdotmapset(d,Ad);
    if (iter==1)
      old_dAd=dAd;
    end
    last_failed=last_failed+1;
    if (abs(dAd)>dad_thresh*abs(old_dAd))
      mdisp(['Damn it.  Failed iteration with ' num2str(dAd) '.  Trying again.']);
      if (last_failed>3)
        mdisp('failed too many times in a row in a row.  unhappily, I am now dead.');
        error('So long, and thanks for all the memory fragmentations.');
      end
     
    else

      last_failed=0;
      alpha=rMr/dAd;
      x=add_mapset(x,d,alpha);
      r=add_mapset(r,Ad,-alpha);
      clear Ad;
      if do_precon,
        Mr=feval(fun,r,precon);
      else
        Mr=r;
      end
      
      rpp=mapsetdotmapset(r,Mr);
      beta=rpp/rMr;
      rMr=rpp;
      if (rMr<best_rMr)
        %best=x;
        best_rMr=rMr;
      end
      d=add_mapset(Mr,d,beta);
      clear Mr;
      iter=iter+1;
      bb=now;

      if (profile)
        fprintf(logid,'%4d %3d %8.2f %8.2f %8.2f %8.2f ',iter,myid,node_time,86400*(cc-aa),86400*(bb-aa),mpi_time);
        for jj=1:length(tods),
          fprintf(logid,'%12.2f %6.2f %6.2f %6.2f %6.2f',tod_times(jj,1),sum(tod_times(jj,2:4)),tod_times(jj,2),tod_times(jj,3),tod_times(jj,4));
        end
        fprintf(logid,'\n');
        fflush(logid);
      end
      
      if (myid==1)
        %disp([iter rMr dAd beta]);

        printf('%4d %14.5g %14.5g %14.5g %9.2f %9.2f %4d %9.2f %9.2f\n',iter,rMr,dAd,beta,86400*(bb-aa),86400*(cc-aa),bad_id,slowest_time,total_time/nproc);
        do_print=false;
        if length(save_interval)==1
          if save_interval>0,
            if rem(iter,save_interval)==0,
              do_print=true;
            end
          end
        else
          if sum(iter==save_interval)>0
            do_print=true;
          end
        end
        if (do_print)
          if isfield(x,'skymap')
            octave2skymap(x.skymap);
            write_map(x.skymap.mapptr,[save_tag num2str(iter)]);
            %write_simple_map_c(x.skymap.mapptr,[save_tag num2str(iter) '.map']);          
          end
        end
      end
      old_dAd=dAd;
      
    end
    
    
    %toc
    
end
if (profile)
  fclose(logid);
end
