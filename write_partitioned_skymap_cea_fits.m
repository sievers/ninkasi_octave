function[value]=write_partitioned_skymap_cea_fits(map,fname,comments)

assert(isstruct(map));
assert(isfield(map,'partition'));
assert(isfield(map,'map'));


rapix=map.partition.global_rapix;
decpix=map.partition.global_decpix;
radelt=map.partition.global_radelt;
decdelt=map.partition.global_decdelt;
pv=map.partition.global_pv;
nx=map.partition.global_nx;
ny=map.partition.global_ny;


keys={};
vals={};

%[keys,vals]=set_keyval_val('SIMPLE','T',keys,vals);
[keys,vals]=set_keyval_val('SIMPLE',true,keys,vals);
[keys,vals]=set_keyval_val('BITPIX',-64,keys,vals);
[keys,vals]=set_keyval_val('NAXIS',2,keys,vals);
[keys,vals]=set_keyval_val('NAXIS1',nx,keys,vals);
[keys,vals]=set_keyval_val('NAXIS2',ny,keys,vals);
%[keys,vals]=set_keyval_val('EXTEND','T',keys,vals);
[keys,vals]=set_keyval_val('EXTEND',true,keys,vals);

%[keys,vals]=set_keyval_val('CTYPE1','''RA---CEA''',keys,vals);
%[keys,vals]=set_keyval_val('CTYPE2','''DEC--CEA''',keys,vals);
%[keys,vals]=set_keyval_val('CUNIT1','''        ''',keys,vals);
%[keys,vals]=set_keyval_val('CUNIT2','''        ''',keys,vals);

[keys,vals]=set_keyval_val('CTYPE1','RA---CEA',keys,vals);
[keys,vals]=set_keyval_val('CTYPE2','DEC--CEA',keys,vals);
[keys,vals]=set_keyval_val('CUNIT1','        ',keys,vals);
[keys,vals]=set_keyval_val('CUNIT2','        ',keys,vals);

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


%if ~is_map_polarized(map.skymap)
mm=map.map; 
if ndims(map.map)==2 %unpolarized map
  write_fits_cell_partitioned(fname,map.map,keys,vals);
else
  for j=1:size(map.map,1)
    ii=max(strfind(fname,'.fits'));
    mytag=['_' get_map_poltag(map.mapptr,j)];
    if ~isempty(ii)
      fname_use=[fname(1:ii-1) mytag fname(ii:end)];
    else
      fname_use=[fname mytag '.fits'];
    end
    write_fits_cell_partitioned(fname_use,squeeze(map.map(j,:,:)),keys,vals);
  end
end

      

