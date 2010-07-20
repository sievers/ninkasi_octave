function[seasons]=guess_tod_season(tod_names)
ctimes=get_tod_ctimes_from_names(tod_names);

starting_ctimes=[0 1217478623 1240359480 1269310955];

seasons=zeros(size(ctimes));
for j=1:length(starting_ctimes),
  ind=ctimes>starting_ctimes(j);
  seasons(ind)=j;
end

return


%+------------+------------+
%| min(ctime) | max(ctime) |
%+------------+------------+
%| 1217478623 | 1232614872 |
%+------------+------------+
%
%Season 3:
%
%+------------+------------+
%| min(ctime) | max(ctime) |
%+------------+------------+
%| 1240359480 | 1263904794 |
%+------------+------------+
%
%Season 4:
%
%+------------+------------+
%| min(ctime) | max(ctime) |
%+------------+------------+
%| 1269310955 | ... |
%+------------+------------+



