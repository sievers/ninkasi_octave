function[value]=write_skymap_radec_fits(map,fname,comments)
[ramin,ramax,decmin,decmax,pixsize]=get_skymap_rect_params_c(map);
ramin=ramin*180/pi;
decmin=decmin*180/pi;
pixsize=pixsize*180/pi;

mm=skymap2octave(map);

keys={};
vals={};

if length(size(mm)==3)
  n1=size(mm,2);
  n2=size(mm,3);
else
  n1=size(mm,1);
  n2=size(mm2);
end

[keys,vals]=set_keyval_val('SIMPLE','T',keys,vals);
[keys,vals]=set_keyval_val('BITPIX',-64,keys,vals);
[keys,vals]=set_keyval_val('NAXIS',2,keys,vals);
[keys,vals]=set_keyval_val('NAXIS1',n1,keys,vals);
[keys,vals]=set_keyval_val('NAXIS2',n2,keys,vals);
[keys,vals]=set_keyval_val('EXTEND','T',keys,vals);

[keys,vals]=set_keyval_val('CTYPE1','''RA---CAR''',keys,vals);
[keys,vals]=set_keyval_val('CTYPE2','''DEC--CAR''',keys,vals);
[keys,vals]=set_keyval_val('CUNIT1','''        ''',keys,vals);
[keys,vals]=set_keyval_val('CUNIT2','''        ''',keys,vals);

[keys,vals]=set_keyval_val('CRVAL1',0,keys,vals);
[keys,vals]=set_keyval_val('CRVAL2',0,keys,vals);
[keys,vals]=set_keyval_val('CDELT1',pixsize,keys,vals);
[keys,vals]=set_keyval_val('CDELT2',pixsize,keys,vals);
[keys,vals]=set_keyval_val('CRPIX1',-ramin/pixsize,keys,vals);
[keys,vals]=set_keyval_val('CRPIX2',-decmin/pixsize,keys,vals);
%[keys,vals]=set_keyval_val('PV2_1',pv,keys,vals);
[keys,vals]=set_keyval_val('EQUINOX',2000,keys,vals);

if exist('comments')
  if ischar(comments)
    comments={comments};
  end
  for j=1:length(comments),
    keys(end+1)='comment';
    vals(end+1)=comments{j};
  end
end

if ~is_map_polarized(map)
  write_fits_cell(fname,mm,keys,vals);
else
  for j=1:size(mm,1)
    ii=max(strfind(fname,'.fits'));
    mytag=['_' get_map_poltag(map,j)];
    if ~isempty(ii)
      fname_use=[fname(1:ii-1) mytag fname(ii:end)];
    else
      fname_use=[fname mytag '.fits'];
    end
    write_fits_cell(fname_use,squeeze(mm(j,:,:)),keys,vals);
  end
end


