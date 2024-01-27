function fitting_data
clear all; close all; clc;
global inputs IC u_data model
plot_skip = 5;
patient_start = 1;
patient_end = 142;
% models = {'full_model',...
%     'full_model_noA2',...
%     'FKPP_model',...
%     'FKPP_model_wA'};
models = {'full_model'};
for im = 1:length(models)
    % Setup lower bounds and upper bounds
    if strcmp(models{im},'full_model')
        model = 'full_model';
        % D1          % D2          % r1         % r2
        LB(1) = 0;    LB(2) = 0;    LB(3) = 0;   LB(4) = 0;
        UB(1) = 0.2;  UB(2) = 0.2;  UB(3) = 15;  UB(4) = 15;
        % K1          % K2          % p          % A2
        LB(5) = 0;    LB(6) = 0;    LB(7) = 0;   LB(8) = 0;
        UB(5) = 1;    UB(6) = 1;    UB(7) = 1;   UB(8) = 3;
    elseif strcmp(models{im},'full_model_noA2')
        model = 'full_model_noA2';
        % D1          % D2          % r1         % r2
        LB(1) = 0;    LB(2) = 0;    LB(3) = 0;   LB(4) = 0;
        UB(1) = 0.2;  UB(2) = 0.2;  UB(3) = 15;  UB(4) = 15;
        % K1          % K2          % p
        LB(5) = 0;    LB(6) = 0;    LB(7) = 0;
        UB(5) = 1;    UB(6) = 1;    UB(7) = 1;
    elseif strcmp(models{im},'FKPP_model')
        model = 'FKPP_model';
        % D1          % r1          % K1
        LB(1) = 0;    LB(2) = 0;    LB(3) = 0;
        UB(1) = 0.2;  UB(2) = 15;   UB(3) = 1;
    elseif strcmp(models{im},'FKPP_model_wA')
        model = 'FKPP_model_wA';
        % D1          % r1          % K1        % A1
        LB(1) = 0;    LB(2) = 0;    LB(3) = 0;  LB(4) = 0;
        UB(1) = 0.2;  UB(2) = 15;   UB(3) = 1;  UB(4) = 3;
    else
        disp("Enter wrong model name")
        break;
    end
    
    % Make save directory
    save_dir = strcat('parameters/',model);
    if ~exist(save_dir, 'dir')
       mkdir(save_dir)
    end    
    
    for i=patient_start:patient_end
        % Save file
        save_file = strcat(save_dir,'/params_cell');
        save_file = strcat(save_file,num2str(i),'.mat');
        % if ~exist(save_file,'file')
            
            % Load data file
            data_file = strcat('data/densities_cell',num2str(i));
            data_file = strcat(data_file,'.mat');
            load(data_file)
        
            IC = IC;
            inputs.t = t;
            % x cannot be 0 for spheroid data
            inputs.x = x+x(2)-x(1);
            u_data = densities;
            
            % Setup for DIRECT
            Problem.f = 'objective';
            bounds = [LB', UB'];
        
            options.maxevals  = 1e4;
            options.maxits    = 2e2;
            options.maxdeep   = 1e4;
            options.showits   = 1;
            options.tol       = 0.01;
            
            % Call DIRECT
            [fval_direct,params_direct] = Direct(Problem,bounds,options);
            
            % Run fmincon
            params0 = params_direct;
            modelfun = @(params) objective(params);
            opts = optimset('MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);
            [params_fmin,fval_fmin] = fmincon(modelfun,params0,[],[],[],[],LB,UB,[],opts);
        
            u_sum = run_model(params_fmin);
            
            save(save_file,'fval_direct','params_direct','params_fmin','fval_fmin','u_sum','u_data','IC','x','t','cell_line')
        % end
    end
end

end