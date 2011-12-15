function[chain_small,like_small,multiplicity,derived_small]=compress_chain_chunk(chain,accept,like,derived_params)
%make sure we have the first sample set correctly for bookkeeping

if accept(1)==0,
  accept(1)=1;
end

ind=find(accept==1);
chain_small=chain(ind,:);
if exist('derived_params')
  if ~isempty(derived_params)
    derived_small=derived_params(ind,:);
  else
    derived_small=[];
  end
end

like_small=like(ind);
ind(end+1)=length(like)+1;
multiplicity=diff(ind);
%if accept(end)==1,
%  multiplicity(end+1)=1;
%end

