function[map]=allocate_healpix_map(nside)
map=allocate_ninkasi_skymap(0.1,0,1,0,1);
set_skymap_healpix_ring_c (map,nside);

