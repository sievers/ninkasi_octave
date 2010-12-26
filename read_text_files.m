function[lines]=read_text_files(varargin)
lines={};
for j=1:length(varargin),
  lines=[lines read_text_file(varargin{j})];
end
