


%% CALIBRATION THROUGH NEWTON-GAUSS Non-Linear Least Square Fitting
%-------------------------------------------------------------------------%


% Developer : Gregorio Marchesini 
% Date      : 20 March 2021
% Contact   : gremar@kth.se


% Description :
% ––––––––––––

% The following script inplements a Newton-Method optimization pocedure in
% order to solve a non-linear fitting probelm. In particular the script was
% designed to succesflly accomplish the calibration of a given test
% magnetometer affected da scale factor and offset uncertainties on the
% measurements. The Least Square approach is implemented in order to
% minimize the distance between a given function f(x,b) which models the
% intensity of the geomagnetoic field measured by the magnetometer and the
% real geomagnetic field intensity.

% In the case of this script the model f(x,b) is a function of the
% matrix of observation x [3xn] where n is the number of observations of
% the the magnetic field component [mx my mz]' and the parameters vector 
% that contains scale the scale factors and the offsets for each axis 
% [offx offy offz sx sy sz]

% The Model
% –––––––––

% f(x,b)=((mx-offx)*sx).^2+((my-offy)*sy).^2+((mz-offz)*sz).^2); 

% This model represents the magnetic field intensity obtained from the set
% of measurments after the corretion for offsets and scale factors is
% applied.

% LEAST SQUARE 
% ––––––––––––

% The aim is to minimise the sum of the squares of the residials as in the
% calssical theory of linear square approximation. The minimuzation process
% will output the parameters [offx offy offz sx sy sz]. 

%   [ S  = sim(ri^2)
%   [ ri = H^2-f(x,b)

% AIM

% minimise S -------> dS/db = 0 (find the minimum in the derivative)

% Note that H is the intensity of the magnetic field at the location of the
% calibration. This value MUST BE CHANGED in the code with then intensity
% of the magnetic filed at your location of

% An exaustive intorudution to the problem can be found oon the realtive
% wikipidia page : 
% https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm

% NOTES FROM THE AUTHOR 

% 1) The algorithm implemented is probably alreadu implemented in matlab in
% the optimization toolbox and this script was developed with the idea of
% understanding the method and aquire useful knowledge on the algorithm.
% in case efficency requirements are high, it is adviced to use built-in
% matlab functions 

% 2) It is extremely impostant that the measurment used for the fitting are
% equally distributed around the center, otherwise the function will likely
% fail to obtaine any derired calibartion. This is because the process will
% be stucked in a local minima which is not meaningful at the scope of the
% calibration

% 3) Initail Guess of the parameter vector must be given as input to the
% script. This initial guess can be obtained with a manual calibration at
% first.

% Location of the calibration and IGRF Refernce :
% –––––––––––––––––––––––––––––––––––––––––––––
%
% Location   : Stockholm, KTH University
% Latitude   : 59° 21° 03''
% Longitude  : 18° 04' 23''
% Day        : 18 March 2021  (Thursday)



%% FUNCTION
%------------------------------------------------------------------------%
clear all
close all
load 'mag_batch' % load the random measurments 

% LOAD YOUR OWN OBSERVATIONS 

decimation      = 1/5                            ; % Decimation of the total observed measuremnts 
observation     = mag_batch(1:1/decimation:end,:); % save the measured magnetic field as variable observation after appling the decimation
observation_cor = observation                    ; % corrected observation. It will be updated later

%% Parameters initilization

syms         offx offy offz sx sy sz   ;    % Variables initialization
var_vec    = [offx offy offz sx sy sz] ;    % Vector of the variables 
H          = (51668*10^-9)*10^6        ;    % [muT] microtesla magnetic field intensity at the location of the measurements
it         = 10                        ;    % Number of iterations in the method

%% Model Implementation

f=@(mx,my,mz)((mx/sx-offx)).^2+((my/sy-offy)).^2+((mz/sz-offz)).^2;

% The function is half-symbolic and half-handle. I want to create multiple 
% symbolic functions starting from the observations 

%% Values Initialization

% CHANGE THE FIRST GUESS IF NECESSARY 
% FIRST GUESS  = [offx      offy     offz   sx   sy sz ]

beta_k         = [9 -10 -16  2   2  2]';  % b in the model description
Delta_beta     = ones(6,1);                                       % delta_beta is the parameters step after each iteration  b_k+1 = b_k + delta_b

% NOTE : the offeset are taken from the zero-field measurements that where
% performed as part of the calibration procedure. A isolated cage was used
% to obtaine this initialization parameters. remeber that the offset
% measured in the zero field cage is still affected by the scale factor, so
% I divided the value obtained in the zero-field cage by two as first guess
% considering the scale factor equal to two


%% Multidimensional Function inizialization

f_sub         = f(observation(:,1),observation(:,2),observation(:,3)); % for each vector [mx my mz] I create a symbolic function according to the model
J             = jacobian(f_sub,var_vec);                               % jacobian of the multidimesional function       df1/db1 df1/db2 df1/db3 ..... 
                                                                       %                                                df2/db1 df2/db2 df2/db3 ..... 
                                                                       %                                                df3/db1 df3/db2 df3/db3 .....                              
[n_obs,n_var] = size(J);                                               % check the size of the jacobian

% Implement while loop with a threshold
% –––––––––––––––––––––––––––––––––––––

for ii = 1:it
    
  f_beta_k   =  double(subs(f_sub,var_vec.',beta_k));      % The functions are all evealuated at the guess-value of beta
  DeltaY     =  H^2-f_beta_k                        ;      % This is known factor that comes form the taylor expansion of f(x,b) -----> f(x,b) = f(x,b_k) + J * delta_beta
  J_k        =  double(subs(J,var_vec.',beta_k))    ;      % jacobina in beta_k
  
  % QR decomposition algorithm
  
  % I suggest this website for a complete understanding of the
  % following passages (https://www2.math.uconn.edu/~leykekhman/courses/MATH3795/Lectures/Lecture_8_Linear_least_squares_orthogonal_matrices.pdf)
  
  [m,n]      =  size (J_k)                          ;      % sizof thre matrix
  [Q,R,P]    =  qr(J_k)                             ;      % QR decomposotion
  c          =  Q.'*DeltaY                          ;      % intermidiate parameter
  y          =  R(1:n,1:n) \ c(1:n)                 ;      % intermidiate parameter
  
  % end of the QR algorithm
  
  %develop the QR decomposition in order to enhace stability of the systems
  Delta_beta =  P*y                                  ;     % estimation of the next step   REMEBER :   b_k+1 = b_k + delta_b   
  beta_k     =  beta_k + Delta_beta                  ;

end

% Correction phase
% ––––––––––––––––

observation_cor(:,1) = observation(:,1)*1/beta_k(4)-beta_k(1);
observation_cor(:,2) = observation(:,2)*1/beta_k(5)-beta_k(2);
observation_cor(:,3) = observation(:,3)*1/beta_k(6)-beta_k(3);


%% GARFICAL REPERESENTATION OF THE RESULTS 
%  Check the uniformity of your observations if the method fails 

% Uniform geomagnetic sphere 

[X,Y,Z] = sphere(50);
X       = X*H;
Y       = Y*H;
Z       = Z*H;

figure('Position',[0,0,1200,600]); 
subplot(131)
hold on;
sphere           = mesh(X,Y,Z);
sphere.FaceAlpha = 0.2;
xlabel('X [\muT]')
ylabel('Y [\muT]')
zlabel('Z [\muT]')


scatter3(observation_cor(:,1),observation_cor(:,2),observation_cor(:,3),200,'r','filled')
scatter3(observation(:,1),observation(:,2),observation(:,3),200,'b','filled')
view(30,30)

subplot(132)
percentual_distance_cor=abs(sqrt(sum(observation_cor.^2,2))-H)/H*100;
histogram(percentual_distance_cor,[0:5:300])
xlabel('Percentage Error from real magnetic field intensity (%)')
ylabel('number of observations')
title('After correction')

subplot(133)

percentual_distance_nocor=abs(sqrt(sum(observation.^2,2))-H)/H*100;
histogram(percentual_distance_nocor,[0:5:300])
xlabel('Percentage Error from real magnetic field intensity (%)')
ylabel('number of observations')
title('Before correction')

disp('––––––––––––––––––––––––––––––––––––––––––––––')
disp('––––––––––––––––––––––––––––––––––––––––––––––')

fprintf('offset in x : %g microT\n',beta_k(1))
fprintf('offset in y : %g microT\n',beta_k(2))
fprintf('offset in z : %g microT\n\n',beta_k(3))

fprintf('scale factor in x : %g \n',beta_k(4))
fprintf('scale factor in y : %g \n',beta_k(5))
fprintf('scale factor in z : %g \n',beta_k(6))

disp('––––––––––––––––––––––––––––––––––––––––––––––')
disp('––––––––––––––––––––––––––––––––––––––––––––––')
