function[value]=corrnoise2tod(noise)
assert(size(noise.map,2)==size(noise.vecs,1));

value=noise.map*noise.vecs;
