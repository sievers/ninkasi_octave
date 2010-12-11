function[value]=decimate_vector(vec)

value=0.25*vec(1:2:end-2)+0.5*vec(2:2:end-1)+0.25*vec(3:2:end);
if iseven(length(vec))
  value(end+1)=[vec(end-2)+vec(end-1)*2+vec(end)]/4;
  %mean(vec(end-1:end));
end
