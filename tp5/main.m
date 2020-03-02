clear all;
close all;
geometry = 1; ## 1 - rectangle, 2 - "L"-geometry, 3 - geom_deux_milieux,...
              ## 4 - trois_deux_milieux
size = 10 ## numbre des maillages successifs, maximum 5 pour la geometry1 (!!)
[maillages,h] = maillage(geometry, size);
representation = 0; ## = 1, si on veut les dessigns des vecteurs propres
verification = 0; ## = 1, si on veut calculer les erreures relatives
##[lambda,erreur_norm] = valeurs_propres(maillages, h, representation, verification)
T = 2;
principal_ondes(maillages,h,T)