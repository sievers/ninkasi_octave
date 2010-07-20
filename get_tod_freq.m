function[freq,ar]=get_tod_freq(tod)
if (ischar(tod))
  tod_name=tod;
else
  assert(class(tod)=='int64');
  tod_name=get_tod_name(tod);
end

ind=1:length(tod_name);
fwee=ind(tod_name=='.');


freq='';
ar='';
if isempty(fwee)|max(fwee)==length(tod_name),
  error(['no array extension found in get_tod_freq in tod ' tod_name]);
end

tag=tod_name(max(fwee)+1:end);

if tag(end)=='/'
  tag=tag(1:end-1);
end

if (length(tag)~=3) | (tag(1)~='a') | (tag(2)~='r')
  error(['unrecognized tag ' tag]);
end

if tag(3)=='1'
  freq='145';
end
if tag(3)=='2'
  freq='215';
end
if tag(3)=='3'
  freq='280';
end

if isempty(freq)
  error(['unrecognized array ' tag]);
end

ar=tag;

