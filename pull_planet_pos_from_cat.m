function[mycat]=pull_planet_pos_from_cat(tod,cat)
todname=tod;
if numel(tod)==1
  if class(tod)=='int64'
    todname=get_tod_name(tod);
  end
end
% else
%  todname=tod;
%end

nn=length(cat.todname);
fwee=strsplit(todname,'/',true);
mytail=fwee{end};
mylen=numel(mytail);
for ii=1:nn,
  if strcmp(mytail,cat.todname{ii}),
    ff=fieldnames(cat);
    for j=1:length(ff),
      %to_eval=['mycat.' ff{j} ' = cat.' ff{j} '(' num2str(ii) ');'];
      %disp(to_eval);
      %eval(to_eval);
      eval(['mycat.' ff{j} ' = cat.' ff{j} '(' num2str(ii) ');']);
    end
    return
  end
end