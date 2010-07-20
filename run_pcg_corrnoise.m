function[x,b,r,ax]=run_pcg_corrnoise(tods,x)
%x is a mapset, could be clear

for j=1:length(tods),
    tods(j).data=double(tods(j).data);
end


b=create_initial_mapset(tods,x);
for j=1:length(tods),
    tods(j).data=0;
end

ax=mapset2mapset_corrnoise(tods,x);
r=add_mapset(b,ax,-1);
rr=mapsetdotmapset(r,r);
r0sqr=rr;
iter=1;
tol=1e-4;
d=r;
disp(['r0sqr is ' num2str(r0sqr)]);


while ((rr>r0sqr*tol)&(iter<20)),
    tic
    Ad=mapset2mapset_corrnoise(tods,d);
    dAd=mapsetdotmapset(d,Ad);
    alpha=rr/dAd;
    x=add_mapset(x,d,alpha);
    r=add_mapset(r,Ad,-alpha);
    rpp=mapsetdotmapset(r,r);
    beta=rpp/rr;
    rr=rpp;
    d=add_mapset(r,d,beta);
    iter=iter+1;
    disp([iter rr dAd alpha beta]);
    toc
    
end
