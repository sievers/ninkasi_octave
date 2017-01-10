function[asdf]=write_healpix_ring(fname,map)
mm=skymap2octave(map);

if ~is_map_polarized(map)
  simple_write_healpix(fname,mm,[],'ring');
else
  for j=1:size(mm,1),
    ii=max(strfind(fname,'.fits'));
    %mytag=['_' get_map_poltag(map,j)];
    mytag=['_' get_map_poltag_octave(map,j)];
    if ~isempty(ii)
      fname_use=[fname(1:ii-1) mytag fname(ii:end)];
    else
      fname_use=[fname mytag '.fits'];
    end
    simple_write_healpix(fname_use,squeeze(mm(j,:,:)),[],'ring');
  end
end
