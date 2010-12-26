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
error(['error in set_act_cea - unrecognized array ' ar]);