function[tod]=read_tod_mustang(fname,varargin)

cut_global=get_keyval_default('cut_global',false,varargin{:});
array_expand_fac=get_keyval_default('array_expand_fac',1.0,varargin{:});
scan_expand_fac=get_keyval_default('scan_expand_fac',1.0,varargin{:});

random_reverse=get_keyval_default('random_reverse',false,varargin{:});
random_sign=get_keyval_default('random_sign',false,varargin{:});
scale_fac=get_keyval_default('scale_fac',1.0,varargin{:});

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
  if (random_reverse)
    disp('randomizing direction');
    doflip=round(rand(1));
    if (doflip)
      data=flipud(data);
      errs=flipud(errs);
    end
  end

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

  if array_expand_fac~=1,
    dd=median(dec,2);
    rr=median(ra,2);
    dd=repmat(dd,[1 size(dec,2)]);
    rr=repmat(rr,[1 size(ra,2)]);
    dec=dd+array_expand_fac*(dec-dd);
    ra=rr+array_expand_fac*(ra-rr);
  end

  if scan_expand_fac~=1,
    dd=median(dec,2);
    rr=median(ra,2);
    dd0=mean(dd);rr0=mean(rr);
    dd=repmat(dd,[1 size(dec,2)]);
    rr=repmat(rr,[1 size(ra,2)]);
    ra=ra-rr;
    dec=dec-rr;

    dd=scan_expand_fac*(dd-dd0)+dd0;
    rr=scan_expand_fac*(rr-rr0)+rr0;
    dec=dec+dd;
    ra=ra+rr;
  end

    



  
  dt=median(diff(tt));
  set_tod_dt_c (tod,dt);
  set_tod_ndata_c(tod,size(data,1));
  set_tod_timevec_c(tod,tt);


  set_tod_rowcol_c(tod,rows,cols);
  set_cuts_mustang(tod,errs,rows,cols);
  set_tod_pointing_saved (tod,ra,dec);

  if (random_sign)
    disp('randomizing sign');
    doflip=round(rand(1));
    if (doflip)
      data=-1*data;
    end
  end


  if (scale_fac~=1)
    disp('rescaling data');
    data=data*scale_fac;
  end



  set_tod_data_saved(tod,data);
  set_tod_radec_lims_c (tod);
  tods(j)=tod;
end
tod=int64(tods);
return;
