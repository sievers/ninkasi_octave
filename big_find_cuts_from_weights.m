function[big_cuts]=big_find_cuts_from_weights(wts,rows,cols)

big_cuts=cell(max(rows)+1,max(cols)+1);
for j=1:size(big_cuts,1),
  for k=1:size(big_cuts,2),
    big_cuts(j,k)={[0,2^31-1]};
  end
end

for j=1:size(wts,2),
  crud=find_cuts_from_weights(wts(:,j));
  if (size(crud,1)>1)
    big_cuts(rows(j)+1,cols(j)+1)={crud};
  else
    if (min(crud)~=1)|(max(crud)~=size(wts,1))
      big_cuts(rows(j)+1,cols(j)+1)={crud};
    end
  end
    
end

