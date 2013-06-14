function[vals]=find_all_fft_pads(jmax,varargin)

j=1;
vals=[];
while j<=jmax
  val=find_happy_fft_pad(j,varargin{:});
  vals(end+1)=val;
  j=val+1;
end

