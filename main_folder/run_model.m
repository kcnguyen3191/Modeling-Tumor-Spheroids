function u_sum = run_model(params)
global inputs IC model
%%%%%%%%%%%%%%%%%%% Error/objective function
%%%%%% Inputs 
% params: parameter vector
% inputs: (x,t) values
% IC: initial condition
% model: model name
%%%%%% Output 
% u_sum: the computed simulation
%%%%%%%%%%%%%%%%%%%

% Extract the parameters for the models
if strcmp(model,'full_model')
    D_params = [params(1), params(2)];
    r_params = [params(3), params(4)];
    K_params = [params(5), params(6)];
    weights = [params(7), 1-params(7)];
    A_params = params(8);
    u0 = [weights(1)*IC, weights(2)*IC];
elseif strcmp(model,'full_model_noA2')
    D_params = [params(1), params(2)];
    r_params = [params(3), params(4)];
    K_params = [params(5), params(6)];
    weights = [params(7), 1-params(7)];
    A_params = 0;
    u0 = [weights(1)*IC, weights(2)*IC];
elseif strcmp(model,'FKPP_model')
    D_params = params(1);
    r_params = params(2);
    K_params = params(3);
    A_params = 0;
    u0 = IC;
elseif strcmp(model,'FKPP_model_wA')
    D_params = params(1);
    r_params = params(2);
    K_params = params(3);
    A_params = params(4);
    u0 = IC;
else
    disp("Enter wrong model name")
end


% inputs
x = inputs.x;
t = inputs.t;

% Solve the PDE/ODE. If failed, then assign big values for u
try
    [~, uOut] = ode15s(@(t,u)adv_dif_spherical_RHS(x, D_params, r_params, K_params, A_params, u),t,u0);

    numT = length(t);
    numX = length(x);
    numD = length(D_params);
    if numD == 2
        u_temp = reshape(uOut, [numT, numX, numD]);
        u1 = u_temp(:,:,1);
        u2 = u_temp(:,:,2);
        u_sum = u1 + u2;
    else
        u_sum = uOut;
    end
catch
    u_sum = 1e6*ones(length(t),length(x));
end
end
