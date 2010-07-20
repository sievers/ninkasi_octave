function[x]=run_pcg_corrnoise_precon(tods,x,precon,fun,priorfun,maxiter)
%x is a mapset, could be clear

if exist('precon')
    do_precon=true;
else
    do_precon=false;
end


for j=1:length(tods),
    tods(j).data=double(tods(j).data);
end


b=create_initial_mapset(tods,x);
for j=1:length(tods),
    tods(j).data=0;
end

ax=mapset2mapset_corrnoise(tods,x);
if exist('priorfun')
    ax=feval(priorfun,ax,x);
end

r=add_mapset(b,ax,-1);
if do_precon,
    Mr=feval(fun,r,precon);
else
    Mr=r;
end

d=Mr;


rMr=mapsetdotmapset(r,d);
r0sqr=rMr;
iter=1;
tol=1e-4;
disp(['r0sqr is ' num2str(r0sqr)]);

if ~exist('maxiter')
    maxiter=50;
end


while ((rMr>r0sqr*tol)&(iter<maxiter)),
    tic
    Ad=mapset2mapset_corrnoise(tods,d);
    if exist('priorfun')
        Ad=feval(priorfun,Ad,d);
    end

    dAd=mapsetdotmapset(d,Ad);
    
    alpha=rMr/dAd;
    x=add_mapset(x,d,alpha);
    r=add_mapset(r,Ad,-alpha);
    if do_precon,
        Mr=feval(fun,r,precon);
    else
        Mr=r;
    end
    
    rpp=mapsetdotmapset(r,Mr);
    beta=rpp/rMr;
    rMr=rpp;
    d=add_mapset(Mr,d,beta);
    iter=iter+1;
    disp([iter rMr dAd alpha beta]);
    toc
    
end
