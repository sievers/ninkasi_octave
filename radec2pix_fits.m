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
   disp(['correcting ' num2str(sum(ind)) ' sources.']);
   pix1(ind)=(ra(ind)-360)/params.radelt+params.raoff;
   pix2=params.decoff+sin(dec*pi/180)*180/pi/params.pv/params.decdelt;
 case {'tan'}
  dec0=params.deccent*pi/180;
  ra0=params.racent*pi/180;
  ra=ra*pi/180;
  dec=dec*pi/180;
  cosc=sin(dec0)*sin(dec)+cos(dec0)*cos(dec).*cos(ra-ra0);
  x=cos(dec).*sin(ra-ra0)./cosc;
  
  y=(cos(dec0)*sin(dec)-sin(dec0)*cos(dec).*cos(ra-ra0))./cosc;
  pix1=params.raoff+x/params.radelt*180/pi;
  pix2=params.decoff+y/params.decdelt*180/pi;
  
 otherwise
  error(['Unsupported projection type ' params.fittype ' in radec2pix_fits']);
end

