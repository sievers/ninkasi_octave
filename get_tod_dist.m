function[dists]=get_tod_dist(varargin)
%Find the nearest approach of uncut detectors to a point on the sky for tods.
%If there are three arguments, point in question is interpreted as radians, if there are
%7, then in HH:MM:SS, etc.  Last argument is list of TODs.  If input position is in radians,
%then returned distance is also in radians.  Otherwise in degrees.

if length(varargin)==3
  ra=varargin{1};
  dec=varargin{2};
  tods=varargin{3};
  fac=1;
else
  assert(length(varargin)==7);
  rah=varargin{1};
  ram=varargin{2};
  ras=varargin{3};
  dd=varargin{4};
  dm=varargin{5};
  ds=varargin{6};
  tods=varargin{7};
  ra=(rah+ram/60+ras/3600)*15*pi/180;
  if dd<0
    dec=-(-dd+dm/60+ds/3600)*pi/180;
  else
    dec=(dd+dm/60+ds/3600)*pi/180;
  end
    
  fac=180/pi;
end

dists=zeros(size(tods));
for j=1:length(tods),
  dists(j)=get_tod_dist_c(ra,dec,tods(j));
end
dists=dists*fac;
