function D_matrix = make_Dmat_spherical(x)
%%%%%%%%%%%%%%%%%%% Function to make the second-order differentiation matrix
%%%%%% Inputs 
% x: radial/spatial inputs
%%%%%% Output 
% D_matrix: the computed second-order differentiation matrix
%%%%%%%%%%%%%%%%%%%
delta_r = x(2) - x(1);

% For spheroid data, the radial locations cannot be 0
if x(1) == 0
    x = x+delta_r;
end

x_int = 2:length(x)-1;

r_vec = x;
r_i = r_vec(x_int);
r_ip1 = r_vec(x_int+1);

% Compute the values for the diagonals
upper = r_ip1.^2 ./ r_i.^2;
upper = upper';
mid = -upper - 1;

D_matrix = sparse([x_int x_int x_int]    ,  ...
    [x_int x_int-1 x_int+1],  ...
    [mid ones(size(x_int)) upper]);


r_2 = r_vec(2);
r_1 = r_vec(1);
r_end = r_vec(end);
r_endp1 = r_vec(end) + delta_r;

D_matrix(1,1) = -r_2^2/r_1^2 - 1;
D_matrix(1,2) =  r_2^2/r_1^2 + 1;
D_matrix(length(x),length(x))   = -r_endp1^2 / r_end^2 - 1;
D_matrix(length(x),length(x)-1) =  r_endp1^2 / r_end^2 + 1;
end
