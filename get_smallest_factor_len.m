function[curmin,curval]=get_smallest_factor_len(len,maxcut)
%inputs are len, maxcut.  figure out a length that lines up well with fft's
if ~exist('maxcut'),
    maxcut=floor(0.01*len); %default to no more than 1% of data cut
end

if maxcut<1,
    maxcut=round(maxcut*len);
end

curmin=max(factor(len));
curval=len;
for j=max(1,len-maxcut):len,
    val=max(factor(j));
    if (val<curval),
        curval=val;
        curmin=j;
    end
end

    