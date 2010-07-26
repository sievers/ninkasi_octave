function[tvec]=get_tod_tvec(mytod)
tvec=get_tod_tvec_c(mytod);
if numel(tvec)==1
  tvec=get_tod_dt*(1:get_tod_ndata(mytod))';
end
