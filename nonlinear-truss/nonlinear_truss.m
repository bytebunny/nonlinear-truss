%% NONLINEAR_TRUSS main script to solve CA1.
% The naming convention is adopted from CALFEM manual version 3.4.
% Page and figure numbers refer to Bonet & Wood, 2nd edition. 
% /Rostyslav Skrypnyk

close all
clear variables
format compact
clc

addpath(genpath('~/Documents/MATLAB/calfem-3/')) % Add Calfem routines.

%% Pre-processing
task = 'C'; % Task # (A, B or C) for setting element area and saving files.
angle_factor = 0.25;
if strcmp(task,'A')
    area_factor = 1;
else
    area_factor = 1.05;
    if strcmp(task,'C')
        angle_factor = 0.1025; % tan(alpha) = 1/3;
    end
end
% Choose between hyperelastic material model and 
% elasto-plastic model:
analysis_type = 'elastic'; % 'elastic' or 'plastic'
if strcmp(analysis_type, 'elastic')
    fprintf('ELASTIC analysis was chosen.\n')
elseif strcmp(analysis_type, 'plastic')
    fprintf('ELASTO-PLASTIC analysis was chosen.\n')
else
    error('/// Wrong analysis type value.')
end

save_to_file = false; % Save plot data.

% Geometry of the structure:
alpha = angle_factor * pi; % Truss angle from the task description.
coord = [0     0;
         100 100/tan(alpha);
         200   0]; % CALFEM: Global coordinate matrix, [mm].
dof = [1 2; % dofs of node 1.
       3 4;
       5 6]; % CALFEM: Global topology matrix (list of dofs for each node).
e_dof = [1 1 2 3 4;
         2 3 4 5 6]; % CALFEM: Element topology (connectivity of 
                     % elements and DOFs) matrix.
N_dof = max(max(e_dof(:,2:end)));
N_el = length( e_dof(:,1) ); % Number of elements.
N_el_nodes = 2; % Number of nodes in each element.
N_el_dof = 2; % Number of DOFs per node.
[el_x, el_y] = coordxtr(e_dof, coord, dof, N_el_nodes);

% Load and time stepping:
u_y = -2*norm(coord(1,:)-coord(2,:)); % Displace by twice the initial
                                      % length (see fig. 3.8b on p.90).
dt = 0.0005; % [s]. Many steps are needed when truss areas are different.
total_t = 1; % [s].
N_steps = round(total_t / dt);

% Velocity boundary conditions: velocity is used to conveniently set
% prescribed value by multiplying with time step.
bc = [1 0; % DOF, prescribed value.
      2 0;
      4 u_y / total_t;
      5 0;
      6 0];

dof_free = setdiff( 1:N_dof, bc(:,1), 'stable'  ); % Free DOFs.

% Material and cross-sectional parameters:
params.areas = [0.5; ...% Cross-sectional areas of the bar elements, [mm^2].
                0.5 * area_factor];
params.young = 210.e3; % Young's modulus, [N/mm^2] (see p.89).
params.nu = 0.5; % Poisson's ratio.
params.yield = 25.e3; % Yield stress, [N/mm^2].
params.plast_mod = 1; % Plastic modulus H, [N/mm^2] (smth small compared to Young's modulus).

state_array = struct('plast_strain',repmat({0},1,N_el), ...
                     'harden_param',repmat({0},1,N_el));
state_array_new = state_array;
%state.plast_strain = zeros(N_el,1); % Plastic strain for each element.
%state.harden_param = zeros(N_el,1); % Hardening parameter (accumulated absolute plastic strain) for each element.

% Tolerances:
error_tol = 1.e-6; % Error tolerance for Newton's iteration.
max_iter = 20; % Maximum number of Newton iterations.

%% Processing
% Initialise history and output variables:
u = zeros( numel(dof),1 ); % Global displacement vector.
du = u; % Increment of displacement.
% History containers:
u_hist = zeros( length(dof_free), N_steps ); % Only node 2 can move.
stress_hist = zeros( N_el, N_steps ); % N elements by N steps.
strain_hist = stress_hist;
strain_plast_hist = strain_hist;
vertical_force_hist = 0; % History of DOF 4.

for step=1:N_steps % Time stepping.
    du( bc(:,1) ) = bc(:,2)*dt; % Reason why BC are defined as velocities.
    u_el = extract(e_dof,u); % Element displacements.
    
    %% Equilibrium (Newton's) iteration
    for iter=1:max_iter
        du_el = extract(e_dof,du); % Incremental element displacements. 

        % Initialize sparse matrix triplets I,J,V 
        % for the global stiffness matrix:
        ind = 0; 
        I = zeros(N_el*(N_el_nodes * N_el_dof)^2,1); 
        J = zeros(N_el*(N_el_nodes * N_el_dof)^2,1); 
        V = zeros(N_el*(N_el_nodes * N_el_dof)^2,1);
         
        residual_vector = zeros(N_dof,1);
        
        %% Assemble matrices
        for i=1:N_el % Loop over elements:
            [force,K,...
             stress, strain,...
             state_array_new(i)] = element_routine(el_x(i,:), el_y(i,:), ...
                                           u_el(i,:), du_el(i,:), ...
                                           params, state_array(i), i, ...
                                           analysis_type);           
            % Assemble global stiffness matrix and RHS vector:
            for ii = 1:N_el_nodes * N_el_dof 
                for jj = 1:N_el_nodes * N_el_dof
                    ind = ind+1; 
                    I(ind) = e_dof(i,1+ii); % 1st column in e_dof is element number.
                    J(ind) = e_dof(i,1+jj); 
                    V(ind) = K(ii,jj);
                end
            end
            
            residual_vector(e_dof(i,2:end)) = ...
               residual_vector(e_dof(i,2:end)) + force;
            % Save history:
            stress_hist(i,step) = stress;
            strain_hist(i,step) = strain;
            strain_plast_hist(i,step) = state_array_new(i).plast_strain;
        end
        K_system = sparse(I,J,V);

        %% Solve the system of linear equations:
         delta_du = - K_system(dof_free,dof_free) \ ...
                    residual_vector(dof_free);
        % Update displacement increment:
        du(dof_free) = du(dof_free) + delta_du;
        % Check if (internal - external) forces are zero:
        if norm(residual_vector(dof_free)) <= error_tol
            fprintf('/// Step %d converged in %d iterations.\n',...
                    step, iter)
            break
        end
    end % Newton's iteration.

    %% Update variables:
    u = u + du;
    state_array = state_array_new;
    % Save history:
    u_hist(:,step) = u(dof_free);
    vertical_force_hist = [vertical_force_hist, vertical_force_hist(end) ...
                           - K_system(bc( bc(:,2)<0, 1 ),:) * du]; % Force is minus internal force.

    if iter == max_iter
        error('/// Step %d did NOT converge in %d iterations.',...
              step, max_iter)
    end
end

%% Post-processing
figure(1) % Fig.3.8(b): vertical force vs deflection
initial_length = norm([ el_x(1,2) - el_x(1,1); ... 
                        el_y(1,2) - el_y(1,1) ]);
data_force = [ [0:dt*abs(u_y):abs(u_y)]' / initial_length,...
                vertical_force_hist' / params.young / sum(params.areas) ]; % account for computing for 2 bars instead of 1.
plot(data_force(:,1), data_force(:,2),'o-')  
xlabel('(Y - y) / L, [-]')
ylabel('F / (EA) , [-]')
%xlim([0 2]) % Same limits as in Fig.3.8(b) of Bonet & Wood.
%ylim([-0.15 0.2])
grid on
if save_to_file % Create file for LaTeX.
    f_name = ['../doc/data/force_deflection_', analysis_type, '_', task, '.dat'];
    f_id = fopen(f_name,'w');
    header = '# (Y - y) / L, [-]   F / (EA), [-]';
    
    fprintf(f_id, '%s\n', header);
    fprintf(f_id, '%.4f           %.4f\n', data_force');
    fclose(f_id);
end



figure() % Fig.3.8(c): Kirchhoff stress vs total strain (Constitutive behaviour)
data_stress = [ [ zeros(length(dof_free),1), strain_hist(1,:) ]' ,...
                [ zeros(length(dof_free),1), stress_hist(1,:) *1e-3 ]' ];
plot(data_stress(:,1), data_stress(:,2),'bo-')
xlabel('Strain, [-]')
ylabel('Kirchhoff stress, [kN/mm2]')
xlim([-0.4 0.3]) % Same limits as in Fig.3.8(c).
ylim([-30 30])
grid on
if save_to_file % Create file for LaTeX.
    f_name = ['../doc/data/stress_strain_', analysis_type, '_', task, '.dat'];
    f_id = fopen(f_name,'w');
    header = '# strain, [-]   Kirchhoff stress, [kN/mm2]';
    
    fprintf(f_id, '%s\n', header);
    fprintf(f_id, '%.4f           %.4f\n', data_stress');
    fclose(f_id);
end



if strcmp(analysis_type, 'plastic')
    figure() % Fig.3.8(d): plastic vs total strain
    data_strain = [ [ zeros(length(dof_free),1), strain_hist(1,:) ]' ,...
                    [ zeros(length(dof_free),1), strain_plast_hist(1,:) ]' ];
    plot(data_strain(:,1), data_strain(:,2),'bo-')
    xlabel('Total strain, [-]')
    ylabel('Plastic strain, [-]')
    xlim([-0.4 0.3]) % Same limits as in Fig.3.8(d).
    ylim([-0.27 0.15])
    grid on
    if save_to_file % Create file for LaTeX.
        f_name = ['../doc/data/tot_pl_strain_', task,'.dat'];
        f_id = fopen(f_name,'w');
        header = '# total strain, [-]   plastic strain, [-]';
    
        fprintf(f_id, '%s\n', header);
        fprintf(f_id, '%.4f           %.4f\n', data_strain');
        fclose(f_id);
    end
end



if ~strcmp(task,'A')
    figure() % Upper joing movement
%     plot(coord(1:2,1), coord(1:2,2),'bo--') % Plot left truss.
%     hold on
%     data_motion = [ [coord(2,1) + [ zeros(length(dof_free),1), u_hist(1,:) ],coord(1,1)]', ...
%                     [coord(2,2) + [0:dt*u_y:u_y], coord(1,2)]' ];
    data_motion = [ [ zeros(length(dof_free),1), u_hist(1,:) ]', ...
                    [0:dt*u_y:u_y]' ];
    plot(data_motion(:,1),data_motion(:,2),'bo-')
    if save_to_file % Create file for LaTeX.
        f_name = ['../doc/data/node_displ_', analysis_type, '_', task, '.dat'];
        f_id = fopen(f_name,'w');
        header = '# horiz displ, [mm]   vert displ, [mm]';
    
        fprintf(f_id, '%s\n', header);
        fprintf(f_id, '%.4f           %.4f\n', data_motion');
        fclose(f_id);
    end
end