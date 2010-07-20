
[map,ra,dec]=read_simple_map('temp_map_wdownhill_sim_input.map', 'double');

[map3,ra,dec]=read_simple_map('temp_map_wdownhill_sim_3.map', 'double');
[map10,ra,dec]=read_simple_map('temp_map_wdownhill_sim_10.map', 'double');
[map30,ra,dec]=read_simple_map('temp_map_wdownhill_sim_30.map','double');

mask=map3~=0;



clf;
subplot(2,2,1)
imagesc(ra*180/pi,dec*180/pi,map.*mask,[-500 500]);
set(gca,'ydir','normal','xdir','reverse');
title('Input Sim');
colorbar
subplot(2,2,2);
imagesc(ra*180/pi,dec*180/pi,(map3-map).*mask,[-500 500]);
set(gca,'ydir','normal','xdir','reverse');
title('Error - 3 Large iterations');

subplot(2,2,3);
imagesc(ra*180/pi,dec*180/pi,(map10-map).*mask,[-500 500]);
set(gca,'ydir','normal','xdir','reverse');
title('Error - 10 Large iterations');

subplot(2,2,4);
imagesc(ra*180/pi,dec*180/pi,(map30-map).*mask,[-500 500]);
set(gca,'ydir','normal','xdir','reverse');
title('Error - 30 Large iterations');


print -depsc2 -r300 /cita/d/www/home/sievers/stuff/map60





temp_map_giantprior_wdownhill_10
