function[asdf]=read_all_tod_noises(froot,tods)


froot(end+1)='/';
for j=1:numel(tods),
  tt=get_tod_tags_from_names(get_tod_name(tods(j)));
  read_tod_noise_banded_projvec(tods(j),[froot 'noise_' tt]);
end
