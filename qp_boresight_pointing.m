function[ra,dec,psi]=qp_boresight_pointing(az,el,ct,varargin)
lat=get_keyval_default('lat',[],varargin{:});
lon=get_keyval_default('lon',[],varargin{:});
pitch=get_keyval_default('pitch',0,varargin{:});
roll=get_keyval_default('roll',0,varargin{:});

%if we didn't specify latitude/longitude, see if we set an observatory by name
if isempty(lat) & isempty(lon)
  [lat,lon]=get_obs_location(varargin{:});
  %lat=lat*pi/180;
  %lon=lon*pi/180;
end

if numel(pitch)==1
  pitch(1:numel(az))=pitch;
end

if numel(roll)==1
  roll(1:numel(az))=roll;
end

if numel(lat)==1
  lat(1:numel(az))=lat;
end

if numel(lon)==1
  lon(1:numel(az))=lon;
end


[ra,dec,psi]=qp_boresight_pointing_c(az,el,pitch,roll,lon,lat,ct);


