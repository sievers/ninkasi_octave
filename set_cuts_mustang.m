function[value]=set_cuts_mustang(tod,errs,rows,cols)

crap=big_find_cuts_from_weights(errs,rows,cols);
alloc_tod_cuts_c(tod);

for j=1:size(crap,1),
  for k=1:size(crap,2),
    mycuts=crap{j,k};
    for m=1:size(mycuts,1),
      cuts_extend_c(tod,mycuts(m,1),mycuts(m,2),j-1,k-1);
    end
  end
end


