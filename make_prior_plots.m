function[value]=make_prior_plots()
[map,ra,dec]=read_simple_map('temp_map_pcg_prior_x1_13.map','double');
whos
imagesc(ra*180/pi,dec*180/pi,13000*fliplr(flipud((real(smooth_image(map,3))'))),[-500 500]);set(gca,'ydir','normal','xdir','reverse');colorbar
title('PCG Prior x1, 520 iter');
print -dpng pipeline_images/pcg_prior_x1_520



[map,ra,dec]=read_simple_map('temp_map_pcg_prior_x0_520,map','double');
imagesc(ra*180/pi,dec*180/pi,real(13000*fliplr(flipud((smooth_image(map,3))'))),[-500 500]);set(gca,'ydir','normal','xdir','reverse');colorbar
title('PCG No Prior, 520 iter');
print -dpng pipeline_images/pcg_prior_x0_520


[map,ra,dec]=read_simple_map('temp_map_pcg_prior_x10_13.map','double');
imagesc(ra*180/pi,dec*180/pi,real(13000*fliplr(flipud((smooth_image(map,3))'))),[-500 500]);set(gca,'ydir','normal','xdir','reverse');colorbar
title('PCG Prior x10, 520 iter');
print -dpng pipeline_images/pcg_prior_x10_520


[map,ra,dec]=read_simple_map('temp_map_pcg_prior_x0.3_520,map','double');
imagesc(ra*180/pi,dec*180/pi,real(13000*fliplr(flipud((smooth_image(map,3))'))),[-500 500]);set(gca,'ydir','normal','xdir','reverse');colorbar
title('PCG Prior x0.3, 520 iter');
print -dpng pipeline_images/pcg_prior_x0.3_520.png


[map,ra,dec]=read_simple_map('temp_map_pcg_prior_x10_find5modes_120.map','double');
imagesc(ra*180/pi,dec*180/pi,real(13000*fliplr(flipud((smooth_image(map,3))'))),[-500 500]);set(gca,'ydir','normal','xdir','reverse');colorbar
title('PCG Prior x10, Find 5 modes, 1200 iter');
print -dpng pipeline_images/pcg_prior_x10_5mode_120.png



[map,ra,dec]=read_simple_map('temp_map_pcg_prior_x0_find5modes_120.map','double');
imagesc(ra*180/pi,dec*180/pi,real(13000*fliplr(flipud((smooth_image(map,3))'))),[-500 500]);set(gca,'ydir','normal','xdir','reverse');colorbar
title('PCG Prior x0, Find 5 modes, 120 iter');
print -dpng pipeline_images/pcg_prior_x0_5mode_120.png



