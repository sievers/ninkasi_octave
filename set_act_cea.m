function[value]=set_act_cea(map,ar)
if ~exist('ar')
  ar='ar1';
end
mdisp(ar)
if strcmp(ar,'ar1')
  set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,8554.0,10034.0,0.353814121044493,14055,1714);
  return
end
if strcmp(ar,'ar1_strip_3year')
    set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,10493,9804,0.353814121044493,18381,1284);
    return
end
if strcmp(ar,'ar2')
  set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,8495.0,9909.0,0.353814121044493,13063,1596);
  return
end
if strcmp(ar,'ar2_new')
  set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,9922.0,9750.0,0.353814121044493,14205,1196);
  return
end
if strcmp(ar,'ar1equ')
  set_skymap_cea_predef_c(map,-0.0008334, 0.0008334,7584.0,323.0,1.0,28830,687);
  return
end

if strcmp(ar,'ar1equ_2010')
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.008334, 0.008334,8558.0,265.0,1.0,28672,644);
  return
end

if strcmp(ar,'allequ')  %Hasselfield's preferred header for all equatorial maps
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.00825, 0.00825,7972.0,292.0,1.0,38275,611);
  return
end


if strcmp(ar,'allsouth')  %Hasselfield's preferred header for all equatorial maps
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.0138696776687, 0.0138696776687,10582,9833,0.35381412,18787,1177);
  return
end


error(['error in set_act_cea - unrecognized array ' ar]);
