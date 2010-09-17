function[row,col]=get_darkdets_rowcol(tod,dirroot)
if exist('dirroot')
  darkname=guess_darkdets_name(tod,dirroot);
else
  darkname=guess_darkdets_name(tod);
end

if isempty(darkname)
  if class(tod)=='int64'
    todname=get_tod_name(tod);
  else
    assert(ischar(tod))
    todname=tod;
  end
  warning(['no dark detector file found for file ' todname]);
  row=[];
  col=[];
  return
end
detlist=load(darkname);
row=floor(detlist/32);
col=rem(detlist,32);
