function[which_ar]=get_tod_array(tod)
if ischar(tod)
  todname=tod;
else
  assert(class(tod)=='int64');
  todname=get_tod_name(tod);
end
if todname(end)=='/'
  todname=todname(1:end-1);
end

which_ar='';
if length(todname)<length('a.ar1')
  return;
end
if todname(end-3)~='.'
  return;
end
which_ar=todname(end-2:end);



