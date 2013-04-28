function[found]=read_cuts(tods,dirroot,varargin)
if ~exist('dirroot')
  dirroot='';
end
read_octave=false;
if ~isempty(varargin)
  read_octave=get_keyval_default('read_octave',false,varargin{:});
end


found=true(size(tods));
for j=1:length(tods)

  myroot=get_mydirroot(dirroot,tods(j));
  if isempty(myroot)
    cutname=guess_cuts_name(tods(j));
  else
    cutname=guess_cuts_name(tods(j),myroot);
  end

  %if exist('dirroot')
  %  cutname=guess_cuts_name(tods(j),dirroot);
  %else
  %  cutname=guess_cuts_name(tods(j));
  %end



  if ~isempty(cutname),
    mdisp(['reading cuts ' cutname ' on tod ' get_tod_name(tods(j))]);
    if (read_octave)
      read_cuts_octave(tods(j),cutname);
    else
      read_cuts_c(tods(j),cutname);
    end
  else
    disp(['unable to find cuts on ' get_tod_name(tods(j))]);
    found(j)=false;
  end
  
end
return


function[myroot]=get_mydirroot(dirroot,tod)
myroot='';
if ischar(dirroot)
  myroot=dirroot;
  return
end
if iscell(dirroot)
  myseason=guess_tod_season(get_tod_name(tod));
  myroot=dirroot{myseason};
end
