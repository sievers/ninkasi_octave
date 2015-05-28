function[ra,dec]=parse_radec_text_line(myline)
ispos=find(myline=='+');
isneg=find(myline=='-');

nsign=numel(ispos)+numel(isneg);
if nsign>1
  error(['too many signs in parse_radec_text_line on line ' myline]);
end

if nsign==1,
  ind=[ispos isneg];
  mysign=1;
  if myline(ind)=='-',
    mysign=-1;
  end
  rastr=myline(1:ind-1);
  decstr=myline(ind+1:end);
  tags=safesplit(strtrim(rastr),' :');
  ra=0;
  for j=1:length(tags),
    ra=ra+str2num(strtrim(tags{j}))/60^(j-1);
  end
  ra=ra*15;

  dec=0;
  tags=safesplit(strtrim(decstr),' :');
  for j=1:length(tags),
    dec=dec+str2num(strtrim(tags{j}))/60^(j-1);
  end
  dec=dec*mysign;
else
  tags=safesplit(strtrim(myline),' :');
  nt=length(tags);
  assert(nt<7);
  assert(nt>1);
  
  nra=ceil(nt/2);
  ratags=tags(1:nra);
  dectags=tags(nra+1:end);
  ra=0;
  dec=0;
  for j=1:length(ratags),
    ra=ra+str2num(strtrim(ratags{j}))/60^(j-1);
  end
  ra=ra*15;
  for j=1:length(dectags),
    dec=dec+str2num(strtrim(dectags{j}))/60^(j-1);
  end
end


function[value]=safesplit(line,varargin)
tags=strsplit(strtrim(line),varargin{:});
keep=true(size(tags));
for j=1:length(tags),
  tags{j}=strtrim(tags{j});
  if isempty(tags{j})
    keep(j)=false;
  end
end
value=tags(keep);