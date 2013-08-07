function[value]=debutter_opts(tod,myopts)
if exist('myopts')
  apply_debutter_opts(tod,true,myopts);
else
  apply_debutter_opts(tod,true);
end
