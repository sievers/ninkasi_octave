function[value]=cut_rightleft_cbass(tod,rl,varargin)
if length(tod)>1,
  for j=1:length(tod),
    cut_rightleft_cbass(tod(j),rl,varargin{:});
  end
  return;
end

[left,right]=find_scans_cbass(tod,varargin{:});

if (lower(rl(1))=='r'),
  to_keep=right;
else 
  if (lower(rl(1))=='l')
    to_keep=left;
  else
    error(['unknown direction ' rl ' specified in cut_rightleft_cbass']);
  end
end
if isempty(to_keep),
  warning(['no scans found in tod ' get_tod_name(tod) '.  Are you sure the telescope was scanning?'])
  return;
end

try
  to_cut=invert_regions(tod,to_keep);
catch
  error(['failed invert_regions on tod ' get_tod_name(tod)]);
end


for j=1:size(to_cut,1),
  cut_tod_global_c(tod,to_cut(j,1),to_cut(j,2));
end

function[to_cut]=invert_regions(tod,to_keep);
n=get_tod_ndata(tod);
to_cut=[to_keep(1:end-1,2)+1 to_keep(2:end,1)-1];
if to_keep(1,1)>1,
  to_cut=[2 to_keep(1,1)-1;to_cut];
end

if (to_keep(end,2)<n-1),
  to_cut=[to_cut;to_keep(end,2) n-2];
end

