function[tod_names,myutc]=parse_loic_selected_tods(tod_names,varargin)

froot=get_keyval_default('dir','/project2/r/rbond/actpol/data/season2/merlin',varargin{:});
min_utc=get_keyval_default('min_utc',0,varargin{:});
max_utc=get_keyval_default('max_utc',24,varargin{:});

isok=false(size(tod_names));
myutc=zeros(size(tod_names));
if iscell(tod_names),
  for j=1:length(tod_names),
    %tt=strsplit(tod_names{j},[' ' sprintf('\t')],true);

    tt=strsplit(strtrim(tod_names{j}),[sprintf('\t')],true);
    myutc(j)=str2num(tt{2});
    a=tt{1};
    tmp=[froot '/' a(1:5) '/' a];
    if (strcmp(tt{end},'2'))
      tod_names(j)=tmp;
      isok(j)=true;
    end
  end
end
tod_names=tod_names(isok);
myutc=myutc(isok);
    
if (min_utc<max_utc)
  isok=(myutc>=min_utc)&(myutc<=max_utc);
else
  isok=(myutc>=min_utc)|(myutc<=max_utc);
end
tod_names=tod_names(isok);
myutc=myutc(isok);

