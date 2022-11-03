%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Master SMaRT - Module Robotique                   %
%                  Kit de developpement etudiant                    %
%                                                                   %
% Entree  : (d,w) roues ar echanitllonnees au temps GPS (champ odo) %
% mesures : (GPS.x, GPS.y)                                          %
%                                                                   %
% Auteur : El Badaoui El Najjar Maan 2021                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load data_smart.mat;
% Defintion of parameters needed by the algorithm
% demie-voie = half distance between the 2 rear wheels
e = 0.7325; %Öá³¤
% Radius of the wheels: coefficents computed by an identification filter
Rard =1.985526745599795e-002;Rarg =1.987828776724217e-002; %ÂÖ°ë¾¶

% Adaptation of the measurements of the rear wheels
darg=[0;diff(odo.sarg)];dard=[0;diff(odo.sard)]; %»¡³¤
z2=Rarg*darg;z3=Rard*dard;
%---------------------------------------------------------------------
% Definition of the time vector:
% Get the time and the vector size in the provided dataset 

t= odo.t; %from GPS
n= size(t,1);%lengteh of vector   (t,1)for row ,£¨t,2£©for colo


%------------------------------------------------------------------------
% parameters for saving:
% Initialize here all parameters that you need: 
% the state vector, the matrices P 

% estimates of position

xs=zeros(1,n) % estimation
ys= zeros(1,n);
thetas= zeros(1,n);
Px= zeros(1,n);
Py= zeros(1,n);
Ro_xy= zeros(1,n);
Ptheta= zeros(1,n);
Ro_xtheta= zeros(1,n);
Ro_ytheta= zeros(1,n);
nux= zeros(1,n);%innovation the difference of measurement and predection of measure
nuy= zeros(1,n);

%************************************************************************
% POSITION ESTIMATOR
% localization state : X = [x;y;theta] 
% input of the filter : elementary rotation and translation measured 
%                     by the sensors of rear wheels 
% measurements :      GPS positions
%************************************************************************
% noise of the model: Qalpha

% Define here the matrix of variance-covariance of the error of the model
Qalpha=	[1 0 0;0 1 0;0 0 (0.01*pi)/180];
        
        
%-------------------------------------------------------------------------        
% Definissez ici la matrice d'observation:
C=[1 0 0;0 1 0];

%------------------------------------------------------------------------
%Initialize here the state of your system:

X=[GPS.x(1);GPS.y(1);0]; %starting position

%-------------------------------------------------------------------------
% Initialisez ici la matrice de variance initiale:
P=[GPS.sigma_lon(1)^2 0 0;0 GPS.sigma_lat(1)^2 0;0 0 180/3*pi/180 ];

%-------------------------------------------------------------------------
% erreur_codeur
R2 = 0.0193^2/12; 
% bruit de l'entree
M = [ 0.5 	  0.5;
   	1/(2*e) -1/(2*e)];
Qgamma = R2*M*M';    %input noise(Odo)

% Computation of the Jacobian matrices
% We define the sizes of the matrices of the linear system
A=eye(3,3);   B=zeros(3,2);   
disp('Traitement en cours...');pause(0.1);

%------------------------------------------------------------------
% MAIN LOOP
%------------------------------------------------------------------
for i=1:n, % i = index of odometry measurements 
 % implement the formulas of elementary translation delta (d) and 
 % elementary rotation omega (w)
	d=  ((z2(i)+z3(i))/2); %elementary translation
	w= ((z2(i)-z3(i))/(2*e));%elementary rotation

   %----------------------------
   %  prediction step 

   % Codez la jacobienne
   A(:,1)=[1;0;0];  A(:,2)=[0;1;0] ;   A(1,3)=-d*sin(X(3)+w/2);
   A(2,3)=d*cos(X(3)+w/2);  A(3,3)=1; 
   B(1,1)= cos(X(3)+w/2);  B(1,2)=-d/2*sin(X(3)+w/2);
   B(2,1)=sin(X(3)+w/2) ;  B(2,2)=(1/2)*d*cos(X(3)+w/2);
   B(3,:)=[0 1];  % derivative for input variable based on f(x,u)
   
   % Codez le modele d evolution
   X(1)=X(1)+ d*cos(X(3)+w/2);
   X(2)=X(2)+ d*sin(X(3)+w/2);
   X(3)=X(3)+ w;
   
   % Codez la matrice de variance-covariance
   P= A*P*A'+Qalpha+B*Qgamma*B';   %we need to consider input error into account
   % it is a position provided in realtime to decrease the latency time
   % Save process
   xs(i)=X(1);   ys(i)=X(2);   thetas(i)=X(3);
   Px(i)=P(1,1); Py(i)=P(2,2); Ro_xy(i)=P(1,2)/sqrt(P(1,1)*P(2,2));
   Ptheta(i)=P(3,3);Ro_xtheta(i)=P(1,3)/sqrt(P(1,1)*P(3,3));Ro_ytheta(i)=P(2,3)/sqrt(P(2,2)*P(3,3));
   
   %----------------------------------------------------------
   % phase d'estimation si on n'est pas dans un masquage
      ro_sx_sy=GPS.sigma_lon(i)*GPS.sigma_lat(i)*GPS.ro(i);
      Qbeta=[GPS.sigma_lon(i)^2,ro_sx_sy;ro_sx_sy,GPS.sigma_lat(i)^2];
       K= P * C' * ( C * P * C' + Qbeta)^-1;;%gain de Kalman
      nu=[GPS.x(i);GPS.y(i)]-[X(1);X(2)];%innovation
      X= X+K*nu;
      P=(eye(3)-K*C)*P*(eye(3)-K*C)+K*Qbeta*K';

   nux(i)=nu(1);nuy(i)=nu(2);
end;

save estim_ekf_GPS_odo.mat t xs ys thetas Px Py Ptheta Ro_xy Ro_xtheta Ro_ytheta nux nuy
disp('EKF termine. Un fichier ".mat" a ete cree');pause(0.1);
hold on; plot(GPS.x, GPS.y, '.b'); plot(xs, ys, 'r')
