function[value]=get_map_poltag_octave(map,ind)
tags={'I','Q','U','QQ','QU','UU'};
polstate=get_map_polstate_c(map);
ii=find(polstate);
value=tags{ii(ind)};
