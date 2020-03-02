function principal_ondes(maillages,h,T, type_MM)
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
c = 2 ;

for h_index = 1:size(h)
  h_i = h(h_index);
  t = cputime;
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
    Mel=matM_elem(S1, S2, S3);
    
  if (type_MM == 0)
    for i=1:3
    I = Numtri(l,i);
      for j=1:3
        J = Numtri(l,j);
        MM(I,J) = MM(I,J) + Mel(i,j);
      end
    end
  end
  
  % On fait l'assemblage des matrices et du second membre
  for i=1:3
    I = Numtri(l,i);
       for j=1:3
          J = Numtri(l,j);
          KK(I,J)=KK(I,J)+Kel(i,j);
       end
   end
end % for l

% Calcul de la CFL
if(type_MM == 0)
lambda = eigs(KK,MM);
lambda_max = lambda(1);
else
k = 1/h_i/h_i;
lambda_max = eigs(k*KK,1);
end
delta_t = sqrt(4/lambda_max);
fprintf('Temps final %6.2f s ; Pas de temps (CFL) %10.6f s\n',T,delta_t);
t_end = cputime - t
end % for maillages

% Nombre de pas de temps
% ----------------------
Nb_temps = round(T/delta_t);

% Debut du chronometre
% ---------------
% A l'aide de cputime

% solution a t=0 (U0(x) = 0 et U1(x) = 0)
% --------------
% UU0 = zeros(Nbpt,1);
% UU1 = ....;

% Boucle sur les pas de temps
% ---------------------------
% Uk-1 = UU0, Uk = UU1 et Uk+1 = UU2
for k = 1:Nb_temps

    tk = k*delta_t;
        
    % Calcul du second membre a l'instant tk
    % --------------------------------------
    % 1. Contribution de la source
    % 2. Contribution des iteres precedents
    % .
    % A COMPLETER
    % .

    % Obtention de Uk+1
    % -----------------
    % .
    % A COMPLETER
    % .
        
    % Visualisation 
    % -------------
    affiche(UU2, Numtri, Coorneu, ['Temps = ', num2str(tk)]);
    axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),...
          max(Coorneu(:,2)),-0.0002 0.0003  -0.00001 0.00015]);
    % Mise a jour des iteres
    % ----------------------
    UU0 = UU1;
    UU1 = UU2;    
end

% Arret du chronometre
% --------------------
% A l'aide de cputime
% temps_calcul = ....

fprintf('Temps de calcul %6.2f s\n',temps_calcul);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

