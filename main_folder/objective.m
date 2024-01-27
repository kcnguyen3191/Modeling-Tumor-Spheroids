function error = objective(params)
global u_data
%%%%%%%%%%%%%%%%%%% Error/objective function
%%%%%% Inputs 
% params: parameter vector
% u_data: data
%%%%%% Output 
% error: sum of squares error
%%%%%%%%%%%%%%%%%%%

% Compute the simulated data
u_sum = run_model(params);

[m_data, n_data] = size(u_data);
[m_sim, n_sim] = size(u_sum);

% Check if the sizes are the same
if m_data == m_sim && n_data == n_sim
    error = sum(sum((u_sum - u_data).^2));
% If not, then assign a big error value
else
    error = 1e9;
end
end