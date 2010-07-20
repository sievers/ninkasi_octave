function[value]=read_text_file(fname)
value={};
fid=fopen(fname,'r');
while 1,
      ll=fgetl(fname);
      if ischar(ll)
            value(end+1)={ll};
      else
        return
      end
end

