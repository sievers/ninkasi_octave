function[big_arr]=mpi_concatenate_cell_strings(arr)

mystr=cell2str(arr);
mystr=mpi_concatenate_string(mystr);
big_arr=str2cell(mystr);
return;

function[mystr]=cell2str(arr)
mystr='';
nl=printf('\n');
for j=1:length(arr),
  mystr=[mystr arr{j} nl];
end

return

function[arr]=str2cell(mystr)
arr={};
nl=printf('\n');
j=0;
while (1)
  j=j+1;
  [tok,mystr]=strtok(mystr,nl);
  if isempty(mystr)
    return;
  end

  arr(j)={tok};
  %disp([ num2str(j) ' .' mystr '.'])
end
