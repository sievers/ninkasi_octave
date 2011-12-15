function[fluxes,weights]=extract_aperture_fluxes(map,params,ra,dec,r1,r2,wt)
xvec=(1:size(map,1))';
yvec=(1:size(map,2));
fluxes=0*ra;
weights=0*ra;
nsrc=numel(ra);
pad=ceil(r2)+1;
padvec=(-pad:pad);
if ~exist('wt')
  wt=[];
end

doweight=~isempty(wt);

for j=1:nsrc,
  [xpix,ypix]=radec2pix_fits(ra(j),dec(j),params);
  xcent=round(xpix);
  ycent=round(ypix);
  xx=xvec(xcent+padvec);
  yy=yvec(ycent+padvec);
  mypatch=map(xcent+padvec,ycent+padvec);
  rmat=repmat((xx-xpix).^2,[1 2*pad+1])+repmat((yy-ypix).^2,[2*pad+1 1]);

  ind1=rmat<=r1^2;
  ind2=(rmat>r1^2)&(rmat<r2^2);
  n1=sum(sum(ind1));
  n2=sum(sum(ind2));
  fluxes(j)=sum(sum(mypatch(ind1)))-sum(sum(mypatch(ind2)))*n1/n2;  
  if doweight,
    mypatch=wt(xcent+padvec,ycent+padvec);
    weights(j)=sum(sum(mypatch(ind1)));
  end
end



