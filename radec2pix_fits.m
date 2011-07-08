function[pix1,pix2]=radec2pix_fits(ra,dec,params)
%convert ra/dec into FITS map coordinates using parameters contained in params.
%params can be set from a file using get_fits_projection_params.m  ra/dec in degrees

if ~isfield(params,'fittype')
  error('Need field fittype set inside params in radec2pix_fits.');
end


switch(params.fittype)
  case {'cea'}
   pix1=ra/params.radelt+params.raoff;
   ind=pix1<0;
   %disp(['correcting ' num2str(sum(ind)) ' sources.']);
   pix1(ind)=(ra(ind)-360)/params.radelt+params.raoff;
   pix2=params.decoff+sin(dec*pi/180)*180/pi/params.pv/params.decdelt;
  otherwise
   error(['Unsupported projection type ' params.fittype ' in radec2pix_fits']);
end

