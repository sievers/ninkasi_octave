function[varargout]=read_many_dirfile_channels(tod_name,varargin)

alphabet=('a'+1)-1;
alphabet=char((1:26)+alphabet-1);


if ischar(tod_name),
  crud=__read_stuff(tod_name,alphabet,varargin{:});
else

  crud=__read_stuff(tod_name{1},alphabet,varargin{:});
  tod_lengths=numel(crud.a);
  for j=2:length(tod_name),
    crap=__read_stuff(tod_name{j},alphabet,varargin{:});
    tod_lengths(end+1)=length(crap.a);
    fn=fieldnames(crap);
    for j=1:length(fn),
      to_eval=['crud.' fn{j} '=[crud.' fn{j} ';crap.' fn{j} '];'];
      eval(to_eval);

    end
  end
end


varargout={};
for j=1:length(fieldnames(crud)),
  to_eval=['varargout(' num2str(j) ')=crud.' alphabet(j) ';'];
  eval(to_eval);
end
if exist('tod_lengths')
  varargout(end+1)=tod_lengths;
end




function[crud]=__read_stuff(tod_name,alphabet,varargin)
lhs=['crud.' alphabet(1)];
for j=2:length(varargin),
  lhs=[lhs  ',crud.' alphabet(j)];
end
to_eval=['[' lhs ']=read_many_dirfile_channels_c(tod_name,varargin{:});'];
eval(to_eval);
