clear all;
close all;
geometry = 1; ## 1 - rectangle, 2 - "L"-geometry, 3 - geom_deux_milieux,...
              ## 4 - trois_deux_milieux
size = 1 ## numbre des maillages successifs, maximum 10 pour la geometry1 (!!)
[maillages,h] = maillage(geometry, size);
type_MM = 0; ##type de matrice de masse qu'on utilise (exacte (type_MM = 0) ou condensee)
representation = 1; ## = 1, si on veut les dessigns des vecteurs propres
verification = 0; ## = 1, si on veut calculer les erreures relatives
##[lambda,erreur_norm] = valeurs_propres(maillages, h, representation, verification)
T = 2;
principal_ondes(maillages,h,T,type_MM, representation)