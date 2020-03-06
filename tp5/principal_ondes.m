function principal_ondes(maillages,h,T, type_MM, representation)
% =====================================================
%
% principal_ondes(nom_maillage,T);
%
% Routine pour la mise en oeuvre des EF P1 Lagrange et
% du schema saute-mouton en temps pour la resolution de
% l'equation des ondes homogene, avec condition aux 
% limites de type Neumann sur le maillage nom_maillage.msh
% Equation des ondes
%
% | d^2u/dt^2 - c^2 Delta u = f,   dans Omega, 0<t<T
% |             du/dn = 0,         sur le bord 0<t<T
% | conditions initiales u(0) = u_0 et du/dt(0) = u_1
%
% =====================================================
fprintf('Eq. des ondes homogene (Elts Finis P1 ; schema saute-mouton)\n');

% Valeurs physiques
% -----------------
% Vitesse de propagation 
c = 2;
for h_index = 1:size(h)
  h_i = h(h_index);
  nom_maillage = maillages{h_index}
  [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=...
  lecture_msh(nom_maillage);

  % declarations
  % ------------
  MM = sparse(Nbpt,Nbpt); % matrice de masse
  KK = sparse(Nbpt,Nbpt); % matrice de rigidite

  % boucle sur les triangles
  % ------------------------
  for l=1:Nbtri
      S1=Coorneu(Numtri(l,1),:);
      S2=Coorneu(Numtri(l,2),:);
      S3=Coorneu(Numtri(l,3),:);
    % calcul des matrices elementaires du triangle l 
    
  % calcul des matrices elementaires (rigidite et masse) 
      Kel=matK_elem(S1, S2, S3, Reftri(l));
       
      if (type_MM == 0) % matrice exacte 
        Mel=matM_elem(S1, S2, S3);
        for i=1:3
        I = Numtri(l,i);
          for j=1:3
            J = Numtri(l,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
            KK(I,J) = KK(I,J)+Kel(i,j);
          end
        end
       else % matrice condensee
    x1 = S1(1); y1 = S1(2);
    x2 = S2(1); y2 = S2(2);
    x3 = S3(1); y3 = S3(2);
    D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)); %l'aire du triangle
    for i=1:3
    I = Numtri(l,i);
    MM(I,I) = MM(I,I) + abs(D)/3;
       for j=1:3
          J = Numtri(l,j);
          KK(I,J)=KK(I,J)+Kel(i,j);
       end
    end
  end %pour if
  
  end % for l
  
  % Calcul de la CFL
  lambda_max = eigs(KK,MM,1);
  delta_t = 2/c/sqrt(lambda_max);
  fprintf('Temps final %6.2f s ; Pas de temps (CFL) %10.6f s\n',T,delta_t)
  
    % Nombre de pas de temps
  % ----------------------
  Nb_temps = round(T/delta_t);

  % Debut du chronometre
  debut_time = cputime;

  % solution a t=0 (U0(x) = 0 et U1(x) = 0)
  F0 = MM*f_gauss(Coorneu(:,1),Coorneu(:,2),0);
  U0 = zeros(Nbpt,1);
  U1 = zeros(Nbpt,1);
  UU0 = U0;
  MUU1 = 0.5*delta_t^2*F0; % MM*UU1 = MM*UU0 - c^2*delta_t^2/2*KK*UU0 + delta_t^2/2*F0 + delta_t *MM* U1, mais UU0 et U1 sont nuls
  UU1 = MM\MUU1;
  % Boucle sur les pas de temps
  % ---------------------------
  % Uk-1 = UU0, Uk = UU1 et Uk+1 = UU2
  E = zeros(Nb_temps,1); % energie discrete
  for k = 1:Nb_temps
    
      tk = k*delta_t;
          
      % Calcul du second membre a l'instant tk
      % --------------------------------------
      % 1. Contribution de la source
      tilde_F = f_gauss(Coorneu(:,1),Coorneu(:,2),tk);

      % Obtention de Uk+1
      % -----------------
      %  UU2 = MM^(-1)*(delta_t^2 * F - c^2 * delta_t^2*KK*UU1 +2*UU1 - UU0), avec F = MM*tilde_F;
      UU2 = (tilde_F - c*c*MM\(KK*UU1))*(delta_t^2) + 2*UU1-UU0;
       
      % energie discrete E(k+1/2)
      E(k) = 0.5*dot(MM*(UU2-UU1),(UU2-UU1))/delta_t/delta_t + ...
      0.5*c*c*dot(KK*UU1,UU2);

      % Visualisation 
      % -------------
      if representation

        affiche(UU2, Numtri, Coorneu, ['Temps = ', num2str(tk)]);
        axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),...
              max(Coorneu(:,2)),-0.0002 0.0003  -0.00001 0.00015]);
       endif  
      % Mise a jour des iteres
      % ----------------------
      UU0 = UU1;
      UU1 = UU2;    
  end
  % Arret du chronometre
  % --------------------
  % A l'aide de cputime
  temps_calcul = cputime - debut_time;
  cptime = cputime;
  
  % evolution d'energie en temps 
  figure
  plot(0:delta_t:4,E)
  grid on;
  xlabel("time, s",'FontSize',18);
  ylabel("E",'FontSize',18);
  title('Energie discrete pour la matrice de masse condensee','FontSize',20)
  
  fprintf('Temps de calcul %6.2f s\n',temps_calcul);
end %pour h


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

