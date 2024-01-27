function A_matrix = make_A_matrix(x)
%%%%%%%%%%%%%%%%%%% Function to make the first-order differentiation matrix
%%%%%% Inputs 
% x: radial/spatial inputs
%%%%%% Output 
% A_matrix: the computed first-order differentiation matrix
%%%%%%%%%%%%%%%%%%%
x_int = 2:length(x)-1;
A_matrix = sparse([x_int x_int x_int]    ,  ...
                      [x_int x_int-1 x_int+1],  ...
       [0*ones(size(x_int)) -1*ones(size(x_int)) 1*ones(size(x_int))]);
   %[-1*ones(size(x_int)) 0*ones(size(x_int)) ones(size(x_int))]);
A_matrix(1,1) = 0;%-1;
A_matrix(1,2) = 0;%1;
A_matrix(length(x),length(x)) = 0;%-1;
A_matrix(length(x),length(x)-1) = 0;%1;
end