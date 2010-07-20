function[new_corrnoise]=tod2corrnoise(data,noise)
assert(size(noise.map,2)==size(noise.vecs,1));
assert(size(data,1)==size(noise.map,1));
assert(size(data,2)==size(noise.vecs,2));
new_corrnoise=noise;
new_corrnoise.map=data*(noise.vecs');

