function[outdir,outroot]=get_outname_from_script(fname,bigset,myset,fac)

fid=fopen(fname,'r');
if fid==-1
  error(['script ' fname ' does not exist in get_outname_from_script']);
end

badval=-12365;

rtstr='';
while ~exist('outroot')
  mystr=fgetl(fid);
  if ischar(mystr)

    if find(mystr=='%')
      mystr=mystr(1:min(find(mystr=='%'))-1);
    end

    %disp(['*' mystr '*'])
    %str1=[rtstr ';' mystr '; rtstr=''''']
    %str2=['rtstr=[rtstr '';'' mystr]']
    %eval(str1,str2);

    if isempty(rtstr)
      rtstr=mystr;
    else
      rtstr=[rtstr ' ; ' mystr];
    end

    disp(rtstr)

    if ~isempty(rtstr)
      clear afdsa
      try
        eval(rtstr)
        rtstr='';
      catch
        disp('combining lines')
      end
    end    
  else
    fclose(fid);
    error(['finished script ' fname ' without hitting outname in get_outname_from_script.']);
  end
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[value]=isbad(a,badval)
if isnumeric(a)
  if numel(a)==1
    if a==badval
      value=true;
      disp('a is bad');
      return
    end
  end
end
disp('a is OK')
value=false;
