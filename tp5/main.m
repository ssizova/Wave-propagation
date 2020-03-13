clear all;
close all;
geometry = 1; ## 1 - rectangle, 2 - geom_deux_milieux,...
              ## 3 - geom_trois_milieux
size = 1 ## nombre des maillages successifs, jusqu'a 10 pour geometry 1
[maillages,h] = maillage(geometry, size);
type_MM = 1; ##type de matrice de masse qu'on utilise (exacte (type_MM = 0) ou condensee)
representation = 1; ## = 1, si on veut les dessigns de la solution
verification = 0; ## = 1, si on veut calculer les erreures relatives
T = 4;
principal_ondes(maillages,h,T,type_MM,representation)