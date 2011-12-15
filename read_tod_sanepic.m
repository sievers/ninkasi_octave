function[tod]=read_tod_sanepic(fname,varargin)
if iscell(fname),
  for j=1:length(fname),
    tod(j)=read_tod_sanepic(fname{j},varargin{:});
  end
  tod=int64(tod);
  return
end

band=get_keyval_default('band','PSW',varargin{:});


[data,ra,dec,mask,timevec,names]=read_sanepic(fname);
ind=false(size(names,1),1);
for j=1:size(names,1),
  if strncmp(names(j,:),band,numel(band))
    ind(j)=true;
  end
end

names=names(ind,:);
data=data(:,ind);
ra=ra(:,ind);
dec=dec(:,ind);
mask=mask(:,ind);
big_rows=names(:,4)-64;
big_cols=0*big_rows;
for j=1:length(big_cols),
  big_cols(j)=str2num(names(j,5:end));
end


tod=allocate_tod_c();
dt=median(diff(timevec));
set_tod_dt_c (tod,dt);
set_tod_ndata_c(tod,size(data,1));
set_tod_timevec_c(tod,timevec);

set_tod_rowcol_c(tod,big_rows,big_cols);

set_cuts_mustang(tod,1e7*mask,big_rows,big_cols);
set_tod_pointing_saved (tod,ra*pi/180,dec*pi/180);
set_tod_data_saved(tod,data);
set_tod_radec_lims_c (tod);


