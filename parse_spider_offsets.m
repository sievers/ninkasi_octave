function[offsets]=parse_spider_offsets(fname)
[a,b]=system(['grep x.r..c.. ' fname]);
ind=find(b==sprintf('\n'));
assert(std(diff(ind))==0);
b(ind)=' ';

b=reshape(b,[ind(1) numel(b)/ind(1)])';
b(:,1)=' ';
b(:,3)=' ';
b(:,6)= ' ';
offsets=str2num(b);

%b(b==sprintf('\n'))=' ';


