function[tods,tod_names,ind]=cull_empty_tods(tods,tod_names,mindet,dofree)
if ~exist('mindet')
   mindet=1;
end
if ~exist('dofree')
  dofree=false;
end



ind=true(size(tods));
assert((length(tods)==length(tod_names))|(length(tod_names)==0));

n=length(tods);
for j=1:n,
    [rows,cols]=get_tod_rowcol(tods(j));
    if length(rows)<mindet,
       ind(j)=false;
       if dofree
         erase_tod_c(tods(j));
       end
    end
end

tods=tods(ind);
if ~isempty(tod_names),
   tod_names=tod_names(ind);
end
