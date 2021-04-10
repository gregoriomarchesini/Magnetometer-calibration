%% MAY THE FORCE BE WITH YOU

% Developer : Gregorio Marchesini
% Date      : 04/04/2021
% mail      : gremar@kth.se


% This script implemets a LEAST SQUARE ELLIPSOIDAL FITTING using
% an adjusted least sqaure estimation method ALS.

% All the passages are taken from the reference [11] in :

% 'Complete Triaxis Magnetometer Calibration in the Magnetic Domain'
% 'Vale ́rie Renaudin, Muhammad Haris Afzal, and Ge ́rard Lachapelle'

% The refernce title is :

% 'Consistent least squares fitting of ellipsoids' 
%  I. Markovsky1, A. Kukush2, S. Van Huffel1

% The procedure is computed step by step from pag 189 of the paper

% NOTE: Some small changes are due to inconsistencies in the paper
% formulatons and were modified accordingly.

%-------------------------------------------------------------------------%
clear all
close all
load 'mag_batch'

decimation        = 1                      ;
n                 = 3                      ;
Mag_obs           = mag_batch(1:end,1:n).' ;  % MATRIX OF THE OBSERVATIONS  3xn 
observation       = Mag_obs                ;
sigma             = .0                     ;  % uT^2 worse case scenario obtained from observation
H                 = 51                     ;  % microT REMEBER TO CHANGE 
[rows,ncol]       = size(Mag_obs)          ;

%% STEP 1 : TENSOR INITIALIZATION 
% candidate solution functions

t0  =  @(x)1                              ;
t1  =  @(x)x                              ;
t2  =  @(x)x.^2-sigma.^2                  ;
t3  =  @(x)x.^3-3*x*sigma.^2              ;
t4  =  @(x)x.^4-6*x.^2*sigma^2+3*sigma.^4 ;

t_vec = {t0,t1,t2,t3,t4}   ;                % The functiuons are stacked in a  cell vector
T     = zeros(5,rows,ncol) ;                % Tensor preallocation

for k=1:5                                   % Note the indices must start from 1 for matlab so I go up until 5 with k. the real range is [0 4].           
   t_sel    = t_vec{k}                ;     % For each page k of the matrix T
   T(k,:,:) = t_sel(observation(:,:)) ;     % you have the whole matrix evaluated as a function
end 

%% STEP 2 : matrix M

SHIFT   = [1:n 0]'              ;
ONE     = ones(n+1,1)           ;
M2      = ONE*SHIFT'            ;
M1      = SHIFT*ONE'            ;
M       = [vec_s(M1) vec_s(M2)] ;   % 10x2 matrix (if n=3)
[p,q]   = size(M)               ;

% This matrix represents the indices of an Kroneker product of a vector
% with himself. (look kroneker product)

% the matrix obtained by kroneker product of a vector with himself is
% symmetric. Thus only the upper triangular part can be considered.
% the element M(p,1),M(p,2) indicates the successive indices of the vecor x
% (one observation) once considered as in the upper triangular part

% x (kronker) x = [x1x1  x1x2  x1x3]    if you take p=3 you have 2 2 becuse
%                 [x2x1  x2x2  x2x3]    it is the third element of vec_s of
%                 [x3x1  x3x2  x3x3]    this matrix


%% STEP 3 : Equality operator

% R matrix preallocation
R=zeros(p,p,n);

for ii=1:n
    for pp=1:p
        for qq=pp:p %only the upper triangular part
            R(pp,qq,ii)=double(M(pp,1)==ii)+double(M(pp,2)==ii)+double(M(qq,1)==ii)+double(M(qq,2)==ii);
        end
    end
end
            
%% STEP 4 : Correction factor

ni=zeros(p,p);

for pp=1:p
    for qq=pp:p       % upper triangular
        memory=0;     % it will store the summation
        for ll=1:ncol
            prod=1;   % it will store the product
            
            for ii=1:n
                prod=prod*T(R(pp,qq,ii)+1,ii,ll); %the +1 is required because the idex 0 in T is not supported
            end
            
            memory=memory+prod;
        end
        ni(pp,qq)=memory;
    end
end


%% STEP 5 : Indices set

D1  = 1:(n+1)*n/2        ; % all the indices in vec_s(A)
D2  = ((1:n)+1).*(1:n)/2 ; % indeces of the diagonal elements of A in vec_s(A)
D   = D1                 ; % useful variable 
for g=1:length(D1)         % Inside the for loop the off diagonal elemnts are 
    if any(D1(g)==D2)      % are  marked with a zero and eliminated 
        D(g) = 0         ;
    end
end
D(D==0) = []             ;

%% STEP 6 : PSI matrix formation

% psi matrix preallocation
psi     = zeros(p,p);

for pp=1:p
    for qq=pp:p
        if any(D==pp) && any(D==qq)
             psi(pp,qq)=4*ni(pp,qq) ;
        elseif ~any(D==pp) && ~any(D==qq)
             psi(pp,qq)=ni(pp,qq) ; 
        else
             psi(pp,qq)=2*ni(pp,qq) ; 
        end
    end
end

psi = triu(psi)+triu(psi,1).';       % report the upper triangular part of
                                     % psi in the lower triangular part.
                                     % Now it is symmetric

%% STEP 7 : Smallest Eigenvector 
% The solution is an eigenvalue problem which can be solved strating from
% the matrix psi

[evec,evalue] = eig(psi)           ;
[~,index]     = min(diag(evalue))  ;
beta          = -1 * evec(:,index) ;

% remeber that an eignevector can be multiplied by any scalar and it is still
% an eigenvector fo a given eignevalue

% this is really important because we can have a problem then in the 
% positive difinitness of the matrix Q that is later defined.
% MAKE sure to check if it is the case to change the sign in your case.

%% STEP 8 : Normalization

beta_est = beta/norm(beta)          ; % estimated solution

%% STEP 9 : Reconstruction of the parameters

Q        = vec_s_inv(beta_est(1:6))         ; % Output from the algorithm
u        = beta_est(n*(n+1)/2+1:p-1)        ; % check the main paper to 
k        = beta_est(p)                      ; % to understand the notation

%% Final steps

b           =  -0.5*(Q)^-1*u                 ; % Here the paper has an error. They forgot a minus 
[V,Lambda]  =  eig(Q)                        ;
alpha       = -H^2/(k-b'*V*Lambda*V'*b)      ; % This formula is expanded wrongly in the paper,so it is correctly reported here
A_inv       =  (V*sqrt(alpha*Lambda)*V')     ; % inverse shape matrix
A           =  inv(A_inv)                    ; % shape matrix

%% FINAL CORRECTION

h           =  A_inv*(observation-b)         ; % This is the passage in which you corret the result 

%% Graphics 

figure('Position',[0,0,1200,600]);
[X,Y,Z] = sphere(50)             ;
X       = X*H                    ;
Y       = Y*H                    ;
Z       = Z*H                    ;

hold on                          ;
sphere           = mesh(X,Y,Z)   ;
sphere.FaceAlpha = 0.2           ;

scatter3(observation(1,:),observation(2,:),observation(3,:),200,'b','filled')
scatter3(h(1,:),h(2,:),h(3,:),200,'r','filled');
xlabel('\muT')
ylabel('\muT')
zlabel('\muT')

view(30,30)
axis equal

figure('Position',[0,0,1200,600])
mod_cor  = abs((sqrt((h(1,:).^2+h(2,:).^2+h(3,:).^2))-H))/H*100;
mod_bias = abs((sqrt((observation(1,:).^2+observation(2,:).^2+observation(3,:).^2))-H))/H*100;

subplot(121)
histogram(mod_cor,[0:5:300])
xlabel('percentage relative error on the geomagnetic field magniture (%)')
ylabel('number of observations')
title('After correction')

subplot(122)
histogram(mod_bias,[0:5:300])
xlabel('percentage relative error on the geomagnetic field magniture (%)')
ylabel('number of observations')
title('Before correction')
%% Print
v=diag(A_inv);
disp('––––––––––––––––––––––––––––––––––––––––––––––')
disp('––––––––––––––––––––––––––––––––––––––––––––––')
fprintf('offset in x : %g microT\n',b(1))
fprintf('offset in y : %g microT\n',b(2))
fprintf('offset in z : %g microT\n\n',b(3))

fprintf('scale factor in x : %g \n',v(1))
fprintf('scale factor in y : %g \n',v(2))
fprintf('scale factor in z : %g \n',v(3))
fprintf('\n')
disp('The full correction matrix with non-orthogonslities') 
disp('and soft ferromagnetic interference is')
fprintf('\n')
disp(A_inv)
disp('––––––––––––––––––––––––––––––––––––––––––––––')
disp('––––––––––––––––––––––––––––––––––––––––––––––')


