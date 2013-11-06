
if exist('inited', 'var')
  return
end
parms.dummy=0;

% General parameters
[seed, parms] = take_from_struct(parms, 'seed', 1);


inited = true;