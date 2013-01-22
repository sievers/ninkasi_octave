function[value]=read_text_file(fname)
value={};
if ischar(fname) %if fname is character, read it
  fid=fopen(fname,'r');
  if fid<1
    value={};
    warning(['unable to open file ' fname ' for reading.']);
    return
  end
else
  fid=fname; %passed in a file pointer
end
while 1,
      ll=fgetl(fid);
      if ischar(ll)
            value(end+1)={ll};
      else
        %close if we opened the file
        if ischar(fname)
          fclose(fid);
        end
        return
      end
end

