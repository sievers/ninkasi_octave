function[srccat]=cut_data_near_radec_fast(tod,ra,dec,radius,varargin)
pointing_model=get_keyval_default('pointing_model','actpol',varargin{:});

srccat.ra=ra;
srccat.dec=dec;
srccat.amps=0*ra+1.0;
srccat.oversamp=1;
srccat.inject_sources=true;

nn=100;
srccat.dx=radius/nn;
vec=0*(1:(round(1.1*nn)+1));
vec(1:nn)=1;
srccat.beam=vec';

did_pointing=false;
if strcmp(pointing_model,'actpol')
  precalc_actpol_pointing_exact(tod);
  did_pointing=true;
end

assert(did_pointing);  %if this fails, you need to tell the code how to do your pointing model

allocate_tod_storage(tod);
assign_tod_value(tod,0);

add_srccat2tod(tod,srccat);

dat=get_tod_data(tod);
%dat(dat>0)=1;
aaa=now;
[rows,cols]=get_tod_rowcol(tod);
for j=1:size(dat,2),
  [i1,i2]=cutvec2inds(dat(:,j)>0);
  for k=1:length(i1),
    cuts_extend_c(tod,i1(k)-1,i2(k),rows(j),cols(j));
  end
end
bbb=now;
disp(['took ' num2str(86400*(bbb-aaa)) ' seconds to get cuts.']);
free_tod_storage(tod);
free_tod_pointing_saved(tod);