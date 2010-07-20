function[value]=divide_maps(map,weight)
assert(size(map,1)==size(weight,1));
assert(size(map,2)==size(weight,2));
weight(weight==0)=1;
value=map./weight;
