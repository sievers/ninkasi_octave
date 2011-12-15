function[cal,row,col]=read_calib_factors(fname)
cal=[];
row=[];
col=[];
lines=read_text_file(fname);
if ~isempty(lines)
  eval([lines{1} ';']);
  eval([lines{2} ';']);
  eval([lines{4} ';']);
end

