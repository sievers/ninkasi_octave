function[scanning_samples]=find_nonscanning_ends(az,tol,min_len)
%cut non-scanning ends of a tod.  Defined as az has not moved by a distance tol from the starting/ending positions
if ~exist('tol','var')
  tol=0.1;
end
tol=tol*pi/180;

if ~exist('min_len','var')
  min_len=50; %minimum length, in samples, for an end to be considered non-scanning
end

az=unwrap(az);

i1=min(find(abs(az-az(1))<tol));
i2=max(find(abs(az-az(end))<tol));
disp([min(az) max(az)])

std(az(end-100:end))*180/pi
if i1<min_len,
  i1=1;
end
if (length(az)-i2)<min_len,
  i2=length(az);
end
scanning_samples=[i1 i2];