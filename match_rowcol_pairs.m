function[match]=match_rowcol_pairs(rr,cc,varargin)

matchfile=get_keyval_default('matchfile','/project/r/rbond/sievers/actpol/angles/actpol_ar1_matched_pairs.txt',varargin{:});
cat=load(matchfile);
match=0*rr-1;
for j=1:length(rr),
  myind=find( (rr(j)==cat(:,1))&(cc(j)==cat(:,2)));
  if ~isempty(myind)
    mymatch=find( (rr==cat(myind,3))&(cc==cat(myind,4)));
    if ~isempty(mymatch)
      assert(numel(mymatch)==1);
      match(j)=mymatch;
    end
  else
    myind=find( (rr(j)==cat(:,3))&(cc(j)==cat(:,4)));
    if ~isempty(myind)
      mymatch=find( (rr==cat(myind,1))&(cc==cat(myind,2)));
      if ~isempty(mymatch)
        assert(numel(mymatch)==1);
        match(j)=mymatch;
      end
    end
  end
  
end

