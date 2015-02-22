function[nn]=get_actpol_cutsname(fname,froot,varargin)
fmt=get_keyval_default('date_format','yyyy-mm-dd',varargin{:});


ct=get_tod_ctimes_from_names (fname);



nn=[froot '/' date_from_ctime(ct,fmt) '/' fname];
if exist(nn)
  return;
end

nn=[froot '/' date_from_ctime(ct+86400,fmt) '/' fname];
if exist(nn)
  return;
end

nn=[froot '/' date_from_ctime(ct-86400,fmt) '/' fname];
if exist(nn)
  return;
end
nn='';
warning(['file not found on ' fname '  in ' froot ' with nn ' nn '.']);

