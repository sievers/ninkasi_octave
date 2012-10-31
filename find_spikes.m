function[dat_filt]=find_spikes(tod,varargin)
outer_width=get_keyval_default('width',20,varargin{:});
inner_width=get_keyval_default('core',1,varargin{:});
thresh=get_keyval_default('thresh',8,varargin{:});

nn=get_tod_ndata(tod);
xvec=(0:nn-1)';
xvec(xvec>nn/2)=xvec(xvec>nn/2)-nn;
filt1=exp(-0.5*xvec.^2/outer_width^2);
filt2=exp(-0.5*xvec.^2/inner_width^2);
filt1=filt1/sum(filt1);
filt2=filt2/sum(filt2);
filt=filt2-filt1;
filtft=fft(filt);
dat_org=get_tod_data(tod);
apply_filter_to_tod(tod,filtft);
dat_filt=get_tod_data(tod);
push_tod_data(dat_org,tod);
clear dat_org;


for j=1:size(dat_filt,2),
  spikes=abs(dat_filt(:,j))>thresh*median(abs(dat_filt(:,j)));
  spike_cuts(j)={find(spikes)};
end

[rr,cc]=get_tod_rowcol(tod);

for det=1:length(spike_cuts)
  if length(spike_cuts{det})>0
    vec=spike_cuts{det};
    for j=1:length(vec)
      imin=vec(j)-1;
      imax=imin+2;
      if imax>=nn
        imax=nn-1;
      end
      cuts_extend_c(tod,imin,imax,rr(det),cc(det));
    end
  end
end






