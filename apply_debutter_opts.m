function[value]=apply_debutter_opts(tod,debutter,myopts)
debutter_octave=get_struct_mem(myopts,'debutter_octave');
myexpt=get_struct_mem(myopts,'debutter_expt','act');

if debutter_octave
  debutterworth_octave(tod,debutter,myexpt);
else
  assert(strcmp(myexpt,'act'));
  if debutter
    debutterworth_c(tod);
  else
    rebutterworth_c(tod);
  end
end

