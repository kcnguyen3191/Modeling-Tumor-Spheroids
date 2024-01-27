function du = adv_dif_spherical_RHS(x, D_params, r_params, K_params, A_params, u)
%%%%%%%%%%%%%%%%%%% PDE function
%%%%%% Inputs 
% x: radial/spatial inputs
% D_params: diffusion parameter(s)
% r_params: growth parameter(s)
% K_params: carrying capacity parameter(s)
% A_params: advection parameter
% u: density
%%%%%% Output 
% du: the computed right hand side (RHS)
%%%%%%%%%%%%%%%%%%%

    % Compute the first and second derivative matrices
    D_matrix = make_Dmat_spherical(x);
    A_matrix = make_A_matrix(x);

    dx = x(2) - x(1);
    numX = length(x);
    numD = length(D_params);

    % Compute the RHS for 2-PDE models
    if numD == 2
        % Compute the total population
        u_matrix = reshape(u,[numX, numD]);
        u1 = u_matrix(:,1);
        u2 = u_matrix(:,2);
        u_sum = u1 + u2;

        % Grab the parameters
        D1 = D_params(1);
        D2 = D_params(2);
        r1 = r_params(1);
        r2 = r_params(2);
        K1 = K_params(1);
        K2 = K_params(2);
        A2 = A_params(1);

        % Compute the RHS for each subpop
        du1 = D1*D_matrix*u1/dx^2 + r1*u1.*(1 - u_sum/K1);
        du2 = D2*D_matrix*u2/dx^2 + r2*u2.*(1 - u_sum/K2) - A2*(0.5*A_matrix*u2/dx);
        du = [du1; du2];
    % Compute the RHS for 1-PDE models
    else
        D = D_params;
        r = r_params;
        K = K_params;
        A = A_params;
        u_sum = u;

        % Compute the RHS
        du = D*D_matrix*u/dx^2 - 0.5*A*A_matrix*u/dx + r*u.*(1 - u_sum/K);
    end
end