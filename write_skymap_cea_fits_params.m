function[value]=write_skymap_cea_fits_params(map,params,fname,comments)

%[rapix,decpix,radelt,decdelt,pv]=get_skymap_cea_params_c(map);
radelt=params.radelt;
decdelt=params.decdelt;
rapix=params.raoff;
decpix=params.decoff;
pv=params.pv;

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
if exist('comments')
  if ischar(comments)
    comments={comments};
  end
  for j=1:length(comments),
    keys(end+1)='comment';
    vals(end+1)=comments{j};
  end
end

%this is a test.  this is only a test.  please enjoy this.  I really cannot say how ling this line is.  I do hope it exceeds 80 characters, though!';


write_fits_cell(fname,map,keys,vals);
      

