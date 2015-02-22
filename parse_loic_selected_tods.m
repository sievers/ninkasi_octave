function[tod_names]=parse_loic_selected_tods(tod_names,varargin)

froot=get_keyval_default('dir','/project2/r/rbond/actpol/data/season2/merlin',varargin{:});
isok=false(size(tod_names));
if iscell(tod_names),
  for j=1:length(tod_names),
    tt=strsplit(tod_names{j},[' ' sprintf('\t')],true);
    a=tt{1};
    tmp=[froot '/' a(1:5) '/' a];
    if (strcmp(tt{end},'2'))
      tod_names(j)=tmp;
      isok(j)=true;
    end
  end
end
tod_names=tod_names(isok);

    

