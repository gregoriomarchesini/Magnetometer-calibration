function A=vec_s_inc(v);
% Developer : Gregorio Marchesini 
% Date      : 4 April 2021
% Contact   : gremar@kth.se

% Description
% –––––––––––

% This function converts a vector v of length n(n+1)/2 into a symmetric
% matrix A of dimension n x n
% 

% INPUTS
% ––––––

%   Name     Data Type            Description 

%   v        array                column or raw vector

% OUTPUT
% ––––––––

%   Name     Data Type            Description 

%   A        array                Symmetric matrix A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FUNCTION
m        = -0.5+sqrt(0.25+2*length(v)); % dimesion of the matrix
A        = zeros(m,m);
counter  = 1;

for jj=1:m 
    A(1:jj,jj) = v(counter:counter+jj-1);
    counter    = counter+jj;
end

A        = A+triu(A,1)';
end