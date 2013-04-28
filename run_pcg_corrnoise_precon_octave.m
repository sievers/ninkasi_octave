function[x]=run_pcg_corrnoise_precon_octave(tods,x,b,precon,fun,priorfun,varargin)
%function[x,best]=run_pcg_corrnoise_precon_octave(tods,x,b,precon,fun,priorfun,varargin)

%x is a mapset, could be clear
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;

if (1)
  if numel(varargin)>1,
    clear myopts;
    for j=1:2:numel(varargin)
      eval(['myopts.' varargin{j} ' = varargin{j+1};']);
    end
  else
    myopts=varargin{1};
  end

  maxiter=get_struct_mem(myopts,'maxiter',50);
  tol=get_struct_mem(myopts,'tol',1e-8);
  save_interval=get_struct_mem(myopts,'save_interval',1e-8);
  save_tag=get_struct_mem(myopts,'save_tag','map_');
  dad_thresh=get_struct_mem(myopts,'dad_thresh',1e3);
  write_times=get_struct_mem(myopts,'write_times',false);
  profile=get_struct_mem(myopts,'profile',false);
  cache_iter=get_struct_mem(myopts,'cache_iter',0);
  restart=get_struct_mem(myopts,'restart',false);
  
else
  maxiter=get_keyval_default('maxiter',50,varargin{:});
  tol=get_keyval_default('tol',1e-8,varargin{:});
  save_interval=get_keyval_default('save_interval',1e-8,varargin{:});
  save_tag=get_keyval_default('save_tag','map_',varargin{:});
  dad_thresh=get_keyval_default('dad_thresh',1e3,varargin{:});
  write_times=get_keyval_default('write_times',false,varargin{:});
  profile=get_keyval_default('profile',false,varargin{:});
  cache_iter=get_keyval_default('cache_iter',0,varargin{:});
  restart=get_keyval_default('restart',false,varargin{:});
end



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
        if iscell(tod_name)
          tod_name=tod_name{1};
        end
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


if (cache_iter>0)
  cache_tails={'tic','toc'};
  cache_tag=[save_tag '.cache'];
  cache_raw={[cache_tag '.' cache_tails{1}],[cache_tag '.' cache_tails{2}]};
  
  aa=find(cache_raw{1}=='/');
  fwee=cache_raw{1};
  cache_link={fwee(aa(end)+1:end)};
  aa=find(cache_raw{2}=='/');
  fwee=cache_raw{2};
  cache_link(end+1)={fwee(aa(end)+1:end)};

  if (myid==1)
    cache_raw;
    cache_link;
    system(['mkdir ' cache_raw{1}]);
    system(['mkdir ' cache_raw{2}]);
  end
  mpi_barrier; %this can be needed if the mkdir command is *so* slow that others get to needing it before the mkdir call has completed.
end


%need d,rMr,x,r, old_dAd
just_read=false;
if restart,
  just_read=true;
  %cache_tag=[save_tag '.cache/'];

  try
    iter1=load([cache_raw{1} '/iter'])*load([cache_raw{1} '/completed']);
  catch
    iter1=0;
  end
  try
    iter2=load([cache_raw{2} '/iter'])*load([cache_raw{2} '/completed']);
  catch
    iter2=0;
  end
  if (iter1>iter2)
    cache_tag=cache_raw{1};
    iter=iter1;
  else
    cache_tag=cache_raw{2};
    iter=iter2;
  end
  cache_tag(end+1)='/';
  %iter=load([cache_tag 'iter']);
  mdisp(['restarting on iteration ' num2str(iter) ' using tag ' cache_tag]);

  fid=fopen([cache_tag 'rMr']);
  rMr=fread(fid,1,'double');
  fclose(fid);

  fid=fopen([cache_tag 'old_dAd']);
  old_dAd=fread(fid,1,'double');
  fclose(fid);

  x=read_mapset(tods,[cache_tag 'x']);
  d=read_mapset(tods,[cache_tag 'd']);
  r=read_mapset(tods,[cache_tag 'r']);

  if isfield(b,'skymap')
    d.skymap.mapptr=make_map_copy(b.skymap.mapptr);
    x.skymap.mapptr=d.skymap.mapptr;
    %make_map_copy(b.skymap.mapptr);
    r.skymap.mapptr=d.skymap.mapptr;
    %make_map_copy(b.skymap.mapptr);
    %octave2skymap(d.skymap);
    %octave2skymap(r.skymap);
    %octave2skymap(x.skymap);
  end
  mdisp(['old_dAd and rMr are ' num2str([old_dAd rMr])]);
end



if exist('cache_tag')
  old_cache_tag=cache_tag;
  if old_cache_tag(end)=='/',
    old_cache_tag=old_cache_tag(1:end-1);
  end
end

while ((rMr>r0sqr*tol)&(iter<maxiter)),
  %tic;
  if (profile)  %if profiling, make sure we all start at the same time
    mpi_barrier;  
  end
  if (cache_iter>0)
    if ((rem(iter,cache_iter)==0)||(iter==3))&(just_read==false) %write early on just in case we have a crash.
      mdisp(['Caching map stuff at iteration ' num2str(iter)]);
      ii=1+iseven(round(iter/cache_iter));
      cache_tag=[cache_raw{ii} '/'];
      if (myid==1)
        cache_tag
        %system(['mkdir ' cache_tag ' >& /dev/null']);
        system(['echo 0 > ' cache_tag 'completed']);
      end
      mpi_barrier; %make sure the directory exists before proceeding.
      cache_start=now;
      if (myid==1)
        fid=fopen([cache_tag 'rMr'],'w');
        fwrite(fid,rMr,'double');
        fclose(fid);
        fid=fopen([cache_tag 'old_dAd'],'w');
        fwrite(fid,old_dAd,'double');
        fclose(fid);
        system(['echo ' num2str(iter) ' > ' cache_tag 'iter']);
      end
      save_mapset(x,tods,[cache_tag 'x']);
      save_mapset(d,tods,[cache_tag 'd']);
      save_mapset(r,tods,[cache_tag 'r']);
      mpi_barrier;
      cache_stop=now;
      if (myid==1)
        system(['echo 1 > ' cache_tag 'completed']);
        system(['rm ' old_cache_tag ' >& /dev/null']);
        %system(['ln -s ' cache_raw{ii} ' ' old_cache_tag]);
        %disp(['ln -s ' cache_raw{ii} ' ' old_cache_tag]);
        
        disp(['took ' num2str(86400*(cache_stop-cache_start)) ' to write cache.']);
      end
    end
  end
  just_read=false;

    aa=now;
    %sz1=size(d.skymap.map);
    [Ad,tod_times,node_time,mpi_time]=mapset2mapset_corrnoise_octave(tods,d,varargin{:});
    %mdisp(['Ad^2 is ' num2str(mapsetdotmapset(Ad,Ad))]);
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
          if isfield(x,'srccat')
            write_srccat_mapset(x.srccat,[save_tag num2str(iter) '_srccat.txt']);
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
