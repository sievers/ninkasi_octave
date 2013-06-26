function[map_conv,map2]=fast_2d_convolve(map1,map2,do_deconvolve)
%(de)convolve 2-d maps, padded so that the FFTs are fast.  Note that edge effects will be not necessarily exactly as expected
%map2 is assumed to be the kernel, and so zeros will be inserted into the middle.  If it already matches the padded size of map1, 
%nothing will be done.  If it is already complex, it will be assumed that the (zero-inserted) kernel will have already been FFTed.
persistent pad_lens
if isempty(pad_lens)
  pad_lens=find_all_fft_pads_v2(1e6);
end
if max([length(map1) length(map2)])>max(pad_lens),
  pad_lens=find_all_fft_pads_v2(max([length(map1) length(map2)]));
end

if ~exist('do_deconvolve')
  do_deconvolve=false;
end


size_org=size(map1);

targ_size=[min(pad_lens(pad_lens>=size(map1,1))) min(pad_lens(pad_lens>=size(map1,2)))];
if all(targ_size==size(map1))==0
  %disp('padding map')
  map1(targ_size(1),targ_size(2))=0;
end
map1ft=fft2_r2c(map1);
if (iscomplex(map2)|(all(size(map2)==size(map1ft))));
  assert(all(size(map2)==size(map1ft)));  %make sure sizes agree if we think kernel FFT is coming in
else
  if size(map2,1)<size(map1,1),
    npad=size(map1,1)-size(map2,1);
    nn=floor(size(map2,1)/2);
    map2=[map2(1:nn,:); zeros(npad,size(map2,2));map2(nn+1:end,:)];
  end
  if size(map2,2)<size(map1,2),
    npad=size(map1,2)-size(map2,2);
    nn=floor(size(map2,2)/2);
    map2=[map2(:,1:nn) zeros(size(map2,1),npad) map2(:,nn+1:end)];
  end
  map2=fft2_r2c(map2);
  assert(all(size(map2)==size(map1ft)));
end

if (do_deconvolve)
  map_conv=fft2_c2r(map1ft./map2);
else
  map_conv=fft2_c2r(map1ft.*map2);
end

map_conv=map_conv(1:size_org(1),1:size_org(2));

