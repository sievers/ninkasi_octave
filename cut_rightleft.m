function[value]=cut_rightleft(tods,rl)
if ~exist('rl')
  warning('right/left not specified in cut_rightleft.  returning');
  return
end
assert(ischar(rl));
assert(class(tods)=='int64');
if tolower(rl(1))=='r',
  cut_right=true;
else
  if tolower(rl(1)=='l')
    cut_right=false;
  else
    warning(['unrecognized direction ' rl ' in cut_rightleft.']);
    return;
  end
end
for j=1:length(tods),
  [alt,az]=get_tod_altaz(tods(j));
  turns=find_tod_turnarounds(tods(j));
  for k=1:length(turns)-1,
    if az(turns(k))<az(turns(k+1))
      if cut_right
        cut_tod_global_c(tods(j),turns(k),turns(k+1));
        %disp(['cutting ' num2str([turns(k) turns(k+1)]) ' with az ' num2str([az(turns(k)) az(turns(k+1))])]);
      end
    else
      if ~cut_right
        cut_tod_global_c(tods(j),turns(k),turns(k+1));
        %disp(['cutting ' num2str([turns(k) turns(k+1)]) ' with az ' num2str([az(turns(k)) az(turns(k+1))])]);
      end      
    end
  end
end
