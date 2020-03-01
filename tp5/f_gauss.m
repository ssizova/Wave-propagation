function val = f_gauss(x,y,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_gauss :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f_t(x,y,t)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%       *   t : le temps
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gaussienne centree en 0.1
% ---------
%  t0 = 0.1;
%  val = exp(-8*(t-t0)^2)*exp(-50*((x-3).^2+(y-1).^2));

%Gaussienne centree en -0.2
% ---------
t0 = -0.2;

val = exp(-50*(t-t0)^2)*exp(-50*((x-3).^2+(y-1).^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
