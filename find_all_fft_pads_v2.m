function[vals]=find_all_fft_pads_v2(jmax,pvec)
if ~exist('pvec')
  pvec=[2 3 5];
end

vals=__find_all_fft_pads_v2(jmax,pvec);
vals=sort(vals);

ii=max(find(jmax>vals));
if ii<length(vals)
  vals=vals(1:ii+1);
end

return


function[vals]=__find_all_fft_pads_v2(jmax,pvec)
npow=ceil(log(jmax)/log(pvec(1)));

if numel(pvec)==1,
  vals=ones(1,npow+1);
  for j=1:npow,
    vals(j+1)=vals(j)*pvec;
  end
  return
end

    

ipow=1;
vals=[];
for j=0:npow,
  myvals=ipow*__find_all_fft_pads_v2(ceil(jmax/ipow),pvec(2:end));
  vals=[vals myvals];
  ipow=ipow*pvec(1);
end




