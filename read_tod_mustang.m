function[tod]=read_tod_mustang(fname,varargin)

cut_global=get_keyval_default('cut_global',false,varargin{:});


if iscell(fname)

  for j=1:length(fname)
    %tod(j)=read_tod_mustang(fname{j});
    if j>1
      tod=[tod read_tod_mustang(fname{j},varargin{:})];
    else
      tod=read_tod_mustang(fname{j},varargin{:});
    end
  end
  return
end


[big_rows,big_cols,big_tt,big_errs,big_data,big_ra,big_dec]=convert_mustang_many_scan(fname);

ntod=length(big_rows);
tods=[];
for j=1:ntod,
  tod=allocate_tod_c();
  rows=big_rows{j};
  cols=big_cols{j};
  tt=big_tt{j};
  errs=big_errs{j};
  data=big_data{j};
  ra=big_ra{j};
  dec=big_dec{j};

  ind=true(size(rows));
  for jj=1:length(ind),
    if min(errs(:,jj))>9.9e5
      ind(jj)=false;
    end
  end

  rows=rows(ind);
  cols=cols(ind);
  errs=errs(:,ind);
  data=data(:,ind);
  ra=ra(:,ind);
  dec=dec(:,ind);
  if cut_global,
    ind=not(min(errs,[],2)>9.9e5);
    disp(['cutting ' num2str(sum(ind==false)) ' global samples.']);
    tt=tt(ind);
    errs=errs(ind,:);
    data=data(ind,:);
    dec=dec(ind,:);
    ra=ra(ind,:);
  end



  
  dt=median(diff(tt));
  set_tod_dt_c (tod,dt);
  set_tod_ndata_c(tod,size(data,1));
  set_tod_timevec_c(tod,tt);


  set_tod_rowcol_c(tod,rows,cols);
  set_cuts_mustang(tod,errs,rows,cols);
  set_tod_pointing_saved (tod,ra,dec);
  set_tod_data_saved(tod,data);
  set_tod_radec_lims_c (tod);
  tods(j)=tod;
end
tod=int64(tods);
return;
