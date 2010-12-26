function[ra,dec]=pix2radec_fits(pix1,pix2,params)
%convert pixel coordinate from FITS maps in ra/dec using parameters contained in params.
%params can be set from a file using get_fits_projection_params.m  ra/dec in degrees

if ~isfield(params,'fittype')
  error('Need field fittype set inside params in pix2radec_fits.');
end


switch(params.fittype)
  case {'cea'}
   ra=params.radelt*(pix1-params.raoff);
   dec=asin((pix2-params.decoff)*params.decdelt*params.pv*pi/180)*180/pi;

  otherwise
   error(['Unsupported projection type ' params.fittype ' in radec2pix_fits']);
end

