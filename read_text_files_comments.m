function[lines]=read_text_files_comments(varargin)
lines={};
if iscell(varargin{1})
  ll=varargin{1};
  for j=1:length(ll),
    lines=[lines read_text_file_comments(ll{j},varargin{2:end})];
  end
else
  for j=1:length(varargin),
    lines=[lines read_text_file_comments(varargin{j})];
  end
end