function[newmap,pix1,pix2]=cea2radec(fname,ra,dec)

[map,header,names]=fits_image_read(fname);
whos map
radelt=get_keyword('CDELT1',names,header);
raoff=get_keyword('CRPIX1',names,header);
decdelt=get_keyword('CDELT2',names,header);
decoff=get_keyword('CRPIX2',names,header);
pv=get_keyword('PV2_1',names,header);








%dec=asin(  (decoff-pix2)*decdelt*pv*pi/180)*180/pi;
%ra=radelt*(pix1-raoff)
%disp([dd2dms_string(ra/15) '   ' dd2dms_string(dec)]);

%dec=asin(  (decoff-pix2)*decdelt*pv*pi/180)*180/pi;
%ra=radelt*(pix1-raoff)
pix1=ra/radelt+raoff;
pix2=decoff+sin(dec*pi/180)*180/pi/pv/decdelt;

newmap=interp2(map,pix1,pix2','linear',0);
