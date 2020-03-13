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


% E_h = zeros(size(h),1); % pour evaluation d'energie a chaque pas de h 

for h_index = 1:size(h)
  nom_maillage = maillages{h_index}
  [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=...
  lecture_msh(nom_maillage);
  
  % declarations
  % ------------ 
  MM = sparse(Nbpt,Nbpt); % matrice de masse
  KK_sigma = sparse(Nbpt,Nbpt); % matrice de rigidite

  % boucle sur les triangles: construction de matrices
  % ------------------------
  if(type_MM == 0) % cas de matrice de masse exacte
  for l=1:Nbtri
      S1=Coorneu(Numtri(l,1),:);
      S2=Coorneu(Numtri(l,2),:);
      S3=Coorneu(Numtri(l,3),:);
   
  % calcul de matrices elementaires de rigidite et de masse
    Kel=matK_elem(S1, S2, S3, Reftri(l)); % matrice elementaire de rigidite
    Mel=matM_elem(S1, S2, S3); % matrice elementaire de masse
    
  % On fait l'assemblage des matrices
   for i=1:3
     I = Numtri(l,i);
          for j=1:3
            J = Numtri(l,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
            KK_sigma(I,J) = KK_sigma(I,J)+Kel(i,j);
          end
        end
  end % pour l
  
else   % cas de matrice de masse condensee
  for l=1:Nbtri
      S1=Coorneu(Numtri(l,1),:);
      S2=Coorneu(Numtri(l,2),:);
      S3=Coorneu(Numtri(l,3),:);
   
  % calcul de matrice elementaire de rigidite
    Kel=matK_elem(S1, S2, S3, Reftri(l));
    
  % On fait l'assemblage des matrices
    x1 = S1(1); y1 = S1(2);
    x2 = S2(1); y2 = S2(2);
    x3 = S3(1); y3 = S3(2);
    D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)); % D est, au signe pres, deux fois l'aire du triangle
    for i=1:3
    I = Numtri(l,i);
    MM(I,I) = MM(I,I) + abs(D)/6;
       for j=1:3
          J = Numtri(l,j);
          KK_sigma(I,J)=KK_sigma(I,J)+Kel(i,j);
       end
    end
end % pour l
end % pour if

  % Calcul de la CFL
  lambda_max = eigs(KK_sigma,MM,1);
  delta_t = 2/sqrt(lambda_max);
  fprintf('Temps final %6.2f s ; Pas de temps (CFL) %10.6f s\n',T,delta_t)
  
  % Nombre de pas de temps
  % ----------------------
  Nb_temps = round(T/delta_t);

  E = zeros(Nb_temps,1); % energie discrete
  %U_point = zeros(Nb_temps,1); % la solution en point (6.5,1)
  %point_number = find(((Coorneu(:,1).-6.5).^2 + (Coorneu(:,2).-1.).^2) < 0.002); % nombre du point qui est la plus proche de (6.5,1)
  
  % Debut du chronometre
  debut_time = cputime;

  % solution a t=0 (U0(x) = 0 et U1(x) = 0)
  % ---------------------------------------
  UU0 = zeros(Nbpt,1);
  UU1 = 0.5*delta_t^2*f_gauss(Coorneu(:,1),Coorneu(:,2),0); % M*U1 = M*U0 - delta_t^2/2*K_sigma*U0 + delta_t^2/2*F0 + delta_t*M*U1,...
                                                            % mais U0 et U1 sont nuls et F0 = M*f_gauss(...)
  M_inv_K = MM\KK_sigma;
  
  % Boucle sur les pas de temps
  % ---------------------------
  % Uk-1 = UU0, Uk = UU1 et Uk+1 = UU2
  
  for k = 1:Nb_temps
    
      tk = k*delta_t;
          
      % Calcul du second membre a l'instant tk
      % --------------------------------------
      % Contribution de la source
      tilde_F = f_gauss(Coorneu(:,1),Coorneu(:,2),tk);

      % Obtention de Uk+1
      % -----------------
      
      UU2 = (tilde_F - M_inv_K*UU1)*(delta_t^2) + 2*UU1-UU0; % U2 = M^(-1)*(delta_t^2*F - delta_t^2*K_sigma*U1 +2*U1 - U0), avec F = M*tilde_F;
      
      %maximum_U = max(abs(UU2)) % pour les tests de CFL numerique
       
      % energie discrete E(k+1/2)
      % -----------------
       E(k) = 0.5*dot(MM*(UU2-UU1),(UU2-UU1))/delta_t/delta_t + 0.5*dot(KK_sigma*UU1,UU2);
      
      % la solution en point (6.5,1)
      % -----------------
      % U_point(k) = UU2(point_number); 
      %{
      % Visualisation 
      % -------------
            
       affiche(UU2, Numtri, Coorneu, ['Temps = ', num2str(tk)]);
       axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),...
              max(Coorneu(:,2)),-0.0002 0.0003  -0.00001 0.00015]);
       %}
       
      % Mise a jour des iteres
      % ----------------------
      UU0 = UU1;
      UU1 = UU2;    
  end
  % Arret du chronometre
  % --------------------
  temps_calcul = cputime - debut_time;
  fprintf('Temps de calcul %6.2f s\n',temps_calcul);
  
  affiche(UU2, Numtri, Coorneu, ['Temps = ', num2str(tk)]); % affichage au dernier pas de temps
       axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),...
              max(Coorneu(:,2)),-0.0002 0.0003  -0.00001 0.00015]);
              
  % E_h(h_index) = E(70); % on prend la valeur d'energie apres le moment de stabilisation
  
  
  % evolution d'energie en temps
  % ----------------------------- 
  figure
  plot(0:delta_t:delta_t*(Nb_temps-1),E, 'LineWidth', 3)
  grid on;
  xlabel('Temps, s', 'FontSize', 18);
  ylabel('E', 'FontSize', 18);
  title('Energie discrete pour la matrice de masse condensee', 'FontSize', 20);
   
  %{
  % evolution de solution en point
  % ------------------------------
  figure
  plot(0:delta_t:delta_t*(Nb_temps-1),U_point)
  grid on;
  xlabel('Temps, s', 'FontSize', 18);
  ylabel('U(6.5,1)', 'FontSize', 18);
  title('Evolution de la solution au (6.5,1)', 'FontSize', 20);
 
 %}
  
end %pour h
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

