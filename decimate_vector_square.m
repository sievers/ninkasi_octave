function[value]=decimate_vector_square(vec)
if ~iseven(length(vec))
  vec(end+1)=vec(end);
end
value=0.5*(vec(1:2:end)+vec(2:2:end));
