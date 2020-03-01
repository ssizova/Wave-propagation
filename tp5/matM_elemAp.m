function MelAp = matM_elemAp(S1, S2, S3)

x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));

if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

MelAp = zeros(3,3);

MelAp(1,1) = abs(D)/6;
MelAp(2,2) = MelAp(1,1);
MelAp(3,3) = MelAp(1,1);