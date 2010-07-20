function[tod_names,tag]=read_tod_lists(rootdir,bigset,myset,do_rising,do_setting)


if ~exist('rootdir')
  rootdir='';
end

if isempty(rootdir)
  rootdir='/home/sievers/act/tod_lists/2008_ar1/';
end

if ~exist('do_rising')
  do_rising=true;
end

if ~exist('do_setting')
  do_setting=true;
end

if (~do_rising) & (~do_setting)
  error('c''mon man - you have to do at least one of rising and setting');
end

if ~exist('bigset'),
  bigset=0;
  myset=[0 1 2 3];
end
if isempty(bigset)
  bigset=0;
end


if ~exist('myset')
  warning(['bigset ' num2str(bigset) ' specified, but myset left blank.  Assuming myset 0 through 3.']);
  myset=[0 1 2 3];
end


if isempty(myset)
  myset=[0 1 2 3];
end


if ischar(myset)
  myset=str2num(myset);
end

if isnumeric(bigset)
  bigset=num2str(bigset);
end

tod_names={};
for j=1:length(myset)
  if (do_setting)
    tod_names=[tod_names read_text_file([rootdir '/split' bigset '/todList_' num2str(myset(j)) '_setting.txt'])];
  end
  if (do_rising)
    tod_names=[tod_names read_text_file([rootdir '/split' bigset '/todList_' num2str(myset(j)) '_rising.txt'])];
  end
end

for j=1:length(tod_names),
  tt=tod_names{j};if tt(end)=='/', tod_names(j)={tt(1:end-1)};end;
end

if length(myset)==4,
  tag='all';
else
  tag=['bigset_' num2str(bigset) '_myset'];
  for j=1:length(myset),
    tag=[tag '_' num2str(myset(j))];
  end
end

if (~do_rising)
  tag=[tag '_setting'];
end
if (~do_setting)
  tag=[tag '_rising'];
end
