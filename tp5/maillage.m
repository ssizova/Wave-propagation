function[maillages, h] = maillage(geometry, size)
maillages = cell(1,size);
h = sparse(size,1);
switch geometry
  case 1
    if size == 1
      h = 0.06;
      maillages{1} = ['geomRect_0.06.msh'];
    else
      for i = 1:size
        h(i) = 0.2 - (i - 1)*0.02;
        name = strcat('geomRect_', num2str(h(i)),'.msh');
        maillages{i} = [name];
      endfor
    endif
  case 2
    h = 0.1;
    maillages{1} = ['geom_deux_milieux.msh'];
  case 3
  h = 0.1;
    maillages{1} = ['geom_trois_milieux.msh'];
end