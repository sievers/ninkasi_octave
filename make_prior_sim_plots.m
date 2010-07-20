
[map10,ra,dec]=read_simple_map('temp_map_giantprior_wdownhill_10.map', 'double');
[map30,ra,dec]=read_simple_map('temp_map_giantprior_wdownhill_30.map', 'double');


clf
subplot(2,2,1);
imagesc(ra*180/pi,dec*180/pi,smooth_image(map10,3),[-500 500]);colorbar;
set(gca,'ydir','normal','xdir','reverse');
title('Giant Prior - 10 Large iterations');

subplot(2,2,2);
imagesc(ra*180/pi,dec*180/pi,smooth_image(map30,3),[-500 500]);colorbar;
set(gca,'ydir','normal','xdir','reverse');
title('Giant Prior - 30 Large iterations');





[map10,ra,dec]=read_simple_map('temp_map_smallprior_wdownhill_10.map', 'double');
[map30,ra,dec]=read_simple_map('temp_map_mediumprior_wdownhill_30.map', 'double');

subplot(2,2,3);
imagesc(ra*180/pi,dec*180/pi,smooth_image(map10,3),[-500 500]);colorbar;
set(gca,'ydir','normal','xdir','reverse');
title('No Prior - 10 Large iterations');

subplot(2,2,2);
imagesc(ra*180/pi,dec*180/pi,smooth_image(map30,3),[-500 500]);colorbar;
set(gca,'ydir','normal','xdir','reverse');
title('No Prior - 30 Large iterations');


print -dpng -r300 /cita/d/www/home/sievers/stuff/map60






