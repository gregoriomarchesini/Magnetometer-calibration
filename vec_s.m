function v=vec_s(A);
% Developer : Gregorio Marchesini 
% Date      : 4 April 2021
% Contact   : gremar@kth.se

% Description
% –––––––––––

% This function converts a symmetric matrix A into a column vector 
% containing the upper triangual part of the matrix.

% INPUTS
% ––––––

%   Name     Data Type         Description 

%   A        array(matrix)     symmetric matrix

% OUTPUT
% ––––––––

%   Name     Data Type        Description 

%   v        array(column)    upper triangular part of the matrix A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FUNCTION


[r,c] = size(A);
v     = [];

for i = 1:c;
    v = [v;A(1:i,i)];
end

end
    