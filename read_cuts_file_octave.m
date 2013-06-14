function[cuts]=read_cuts_file_octave(cutsname)
ll=read_text_file(cutsname);
dims=sscanf(ll{1},'%d');
det_cuts=cell(dims(1),dims(2));
for j=1:dims(1),
  for k=1:dims(2)
    det_cuts(j,k)={nan};
  end
end

tags=strsplit(ll{2},'rc:(), ',true);
if ~isempty(tags)
  cuts.global_cuts=cellstr2mat(tags);
else
  cuts.global_cuts=[];
end

for j=3:length(ll),
  tags=strsplit(ll{j},'rc:(), ',true);
  vec=cellstr2mat(tags);
  det_cuts(vec(1)+1,vec(2)+1)={vec(3:end)};

end
cuts.det_cuts=det_cuts;



