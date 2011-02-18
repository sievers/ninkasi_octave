function[poff]=parse_pointing_offsets_type_100(lines)
%return a pointing offset structure from text lines for the keyword-based
%pointing file.  Input should be comment-/whitespace-stripped lines

blocks=extract_blocks(lines);
det_lines=extract_dets(lines);
poff.det_raw=lines2array(det_lines);
poff.pstruct=cell(size(blocks));
for j=1:length(blocks),
  poff.pstruct(j)={block2struct(blocks{j})};
end
return

function[pstruct]=block2struct(lines)
for j=1:length(lines),
  [a,b]=strtok(lines{j});
  bb=sscanf(b,'%f');
  eval(['pstruct.' strtrim(a) ' = bb;']);
end



function[arr]=lines2array(lines);
ll=sscanf(lines{1},'%f')';
nl=length(lines);
arr(nl,length(ll))=0;
arr(1,:)=ll;
for j=2:nl,
  %if this line fails, there are different #'s of entries per line
  arr(j,:)=sscanf(lines{j},'%f')';
end



function[dets]=extract_dets(lines)
for j=1:length(lines),
  if strcmp(lines{j},'begin_offsets');
    ind_start=j;
  end
  if strcmp(lines(j),'end_offsets')
    dets=lines(ind_start+1:j-1);
    return
  end
end
assert(1==0) %never get here, means we didn't have 'end_offsets'
return



function[blocks]=extract_blocks(lines)
nbegin=0;
nend=0;
begin_ind=[];
end_ind=[];
for j=1:length(lines),
  if strcmp(strtrim(lines{j}),'begin_match')
    nbegin=nbegin+1;
    begin_ind(nbegin)=j;
  end
  if strcmp(strtrim(lines{j}),'end_match')
    nend=nend+1;
    end_ind(nend)=j;
  end
end
%disp([ nbegin nend])
assert(length(end_ind)==length(begin_ind));

nblock=length(begin_ind);
blocks=cell(nblock,1);
for j=1:nblock,
  blocks(j)={lines(begin_ind(j)+1:end_ind(j)-1)};
end
return
