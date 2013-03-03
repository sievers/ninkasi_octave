function[value]=write_skymap_cea_fits(map,fname)

if isstruct(map)
  assert(isfield(map,'map'));
    assert(isfield(map,'mapptr'));
    octave2skymap(map);
    map=map.mapptr;
end
assert(class(map)=='int64');

mm=skymap2octave(map);



[rapix,decpix,radelt,decdelt,pv]=get_skymap_cea_params_c(map);
keys={};
vals={};

[keys,vals]=set_keyval_val('SIMPLE','T',keys,vals);
[keys,vals]=set_keyval_val('BITPIX',-64,keys,vals);
[keys,vals]=set_keyval_val('NAXIS',2,keys,vals);
[keys,vals]=set_keyval_val('NAXIS1',size(map,1),keys,vals);
[keys,vals]=set_keyval_val('NAXIS2',size(map,2),keys,vals);
[keys,vals]=set_keyval_val('EXTEND','T',keys,vals);

[keys,vals]=set_keyval_val('CTYPE1','''RA---CEA''',keys,vals);
[keys,vals]=set_keyval_val('CTYPE2','''DEC--CEA''',keys,vals);
[keys,vals]=set_keyval_val('CUNIT1','''        ''',keys,vals);
[keys,vals]=set_keyval_val('CUNIT2','''        ''',keys,vals);

[keys,vals]=set_keyval_val('CRVAL1',0,keys,vals);
[keys,vals]=set_keyval_val('CRVAL2',0,keys,vals);
[keys,vals]=set_keyval_val('CDELT1',radelt,keys,vals);
[keys,vals]=set_keyval_val('CDELT2',decdelt,keys,vals);
[keys,vals]=set_keyval_val('CRPIX1',rapix,keys,vals);
[keys,vals]=set_keyval_val('CRPIX2',decpix,keys,vals);
[keys,vals]=set_keyval_val('PV2_1',pv,keys,vals);
[keys,vals]=set_keyval_val('EQUINOX',2000,keys,vals);
[keys,vals]=set_keyval_val('PC1_1',1.0,keys,vals);
[keys,vals]=set_keyval_val('PC1_2',0.0,keys,vals);
[keys,vals]=set_keyval_val('PC2_1',0.0,keys,vals);
[keys,vals]=set_keyval_val('PC2_2',1.0,keys,vals);

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

      

