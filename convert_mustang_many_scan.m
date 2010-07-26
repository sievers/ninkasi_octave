function[big_rows,big_cols,big_tt,big_errs,big_data,big_ra,big_dec]=convert_mustang_many_scan(fname)
fid=fopen(fname,'r');
[a,b]=read_fits_header(fid);
[a2,b2]=read_fits_hdu(fid);
fclose(fid);

if ~isfield(a2,'SCAN')
  scan=zeros(size(a2.DX));
else
  scan=a2.SCAN;
end

scan_inds=unique(scan);
nscan=length(scan_inds);
for j=nscan:-1:1,
  ind=scan==scan_inds(j);

  npix=length(unique(a2.PIXID(ind)));
  nn=length(a2.DX(ind));
  
  nframe=nn/npix;
  assert(nframe==round(nframe));
  
  tt=reshape(a2.TIME(ind),[nframe,npix]);
  errs=reshape(a2.UFNU(ind),[nframe,npix]);
  rows=reshape(a2.ROW(ind),[nframe,npix]);
  cols=reshape(a2.COL(ind),[nframe,npix]);
  ra=reshape(a2.DX(ind),[nframe,npix]);
  dec=reshape(a2.DY(ind),[nframe,npix]);
  rows=rows(1,:);
  cols=cols(1,:);
  data=reshape(a2.FNU(ind),[nframe,npix]);
  
  tt=tt(:,1);
  tt=rectify_times(tt);

  big_rows(j)={rows};
  big_cols(j)={cols};
  big_tt(j)={tt};
  big_errs(j)={errs};
  big_data(j)={data};
  big_ra(j)={ra};
  big_dec(j)={dec};

end
