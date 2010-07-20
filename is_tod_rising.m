function[isrising]=is_tod_rising(tods)

assert(class(tods)=='int64');

isrising=false(size(tods));
for j=1:length(tods),
  [az,alt,ctime]=get_median_altaz_c(tods(j));
  [ra,dec]=alt_az_ctime2act_ra_dec_c(alt,az,ctime);

  
  %if dec*180/pi>-25
  %  strip_equ='equ';
  %else
  %  strip_equ='south';
  %end

  if (az>pi)
    %riseset='setting';
    isrising(j)=false;
  else
    %riseset='rising';
    isrising(j)=true;
  end
end


