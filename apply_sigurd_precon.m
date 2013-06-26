function[skymap]=apply_sigurd_precon(skymap,precon)
skymap.map=skymap.map.*precon.wtroot;
skymap.map=fast_2d_convolve(skymap.map,precon.impulse);
skymap.map=skymap.map.*precon.wtroot;

