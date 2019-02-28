%% NONLINEAR_TRUSS main script to solve CA1.
% The naming convention is adopted from CALFEM manual version 3.4.
% Page and figure numbers are made to Bonet & Wood, 2nd edition. 
% /Rostyslav Skrypnyk

close all
clear variables
format compact
clc

addpath(genpath('~/Documents/MATLAB/calfem-3/')) % Add Calfem routines.

%% Settings
% Geometry of the structure:
coord = [0     0;
         100 100; % computed via given width and angle of the frame.
         200   0]; % CALFEM: Global coordinate matrix, [mm].
dof = [1 2; % dofs of node 1.
       3 4;
       5 6]; % CALFEM: Global topology matrix (list of dofs for each node).
e_dof = [1 1 2 3 4;
         2 3 4 5 6]; % CALFEM: Element topology (connectivity of 
                     % elements and DOFs) matrix.
N_dof = max(max(e_dof(:,2:end)));
[el_x, el_y] = coordxtr(e_dof, coord, dof, 2);

% Load and time stepping:
u_y = -2*norm(coord(1,:)-coord(2,:)); % Displace by twice the initial
                                      % length (see fig. 3.8b on p.90).
dt = 0.01; % [s].
total_t = 1; % [s].
N_steps = total_t / dt;

% Velocity boundary conditions: velocity is used to conveniently set
% prescribed value by multiplying with time step.
bc = [1 0; % DOF, prescribed value.
      2 0;
      4 u_y / total_t;
      5 0;
      6 0];

dof_free = setdiff( 1:N_dof, bc(:,1), 'stable'  ); % Free DOFs.

% Material and cross-sectional parameters:
params.areas = [1; 1]; % Cross-sectional areas of the bar elements, [mm^2].
params.young = 210.e3; % Young's modulus, [N/mm^2] (see p.89).
params.nu = 0.5; % Poisson's ratio.
params.yield = 25.e3; % Yield stress, [N/mm^2].

% Tolerances:
error_tol = 1.e-6; % Error tolerance for Newton's iteration.
max_iter = 20; % Maximum number of Newton iterations.

%% Solution
% Initialise history and output variables:
u = zeros( numel(dof),1 ); % Global displacement vector.
du = u; % Increment of displacement.
u_hist = zeros( length(dof_free), N_steps ); % Only node 2 can move.
stress_hist = zeros( length(e_dof(:,1)), N_steps ); % N elements by N steps.
vertical_force_hist = 0%zeros( 1, N_steps ); % History of DOF 4.

for step=1:N_steps % Time stepping.
    du( bc(:,1) ) = bc(:,2)*dt; % Reason why BC are defined as velocities.
    u_el = extract(e_dof,u); % Element displacements.
    du_el = extract(e_dof,du); % Incremental element displacements.
    
    %% Equilibrium (Newton's) iteration
    for iter=1:max_iter
% Make sparse later:
        K_system = zeros( numel(dof) ); % Global stiffness
                                                % matrix.
        residual_vector = zeros(numel(dof),1);
        
        %% Assemble matrices
        for i=1:size(e_dof,1) % Loop over elements:
            [force,K,stress] = element_routine(el_x(i,:), el_y(i,:), ...
                                               u_el(i,:), du_el(i,:), ...
                                               params,i);
            % Assemble global stiffness matrix and RHS vector:
            K_system(e_dof(i,2:end),e_dof(i,2:end)) = ...
               K_system(e_dof(i,2:end),e_dof(i,2:end)) + K;
            
            residual_vector(e_dof(i,2:end)) = ...
               residual_vector(e_dof(i,2:end)) + force;
            % Save history:
            stress_hist(i,step) = stress;
        end
        
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
figure() % Vertical force history
initial_length = norm([ el_x(1,2) - el_x(1,1); ... 
                        el_y(1,2) - el_y(1,1) ]);
plot([0:dt*abs(u_y):abs(u_y)] / initial_length, ...
     [ vertical_force_hist / params.young / params.areas(1) / 2 ],'o-') % Reverse sign to match Fig.3.8(b) in Bonet & Wood, 
% and account for computing for 2 bars instead of 1.
xlabel('(Y - y) / L, [-]')
ylabel('F / (EA) , [-]')
xlim([0 2]) % Same limits as in Fig.3.8(b) of Bonet & Wood.
ylim([-0.15 0.2])
grid on

figure() % Constitutive behaviour
plot(0:dt*abs(u_y):abs(u_y), ...
     [ zeros(length(dof_free),1), stress_hist(1,:) *1e-3 ],'bo-')
hold on
plot(0:dt*abs(u_y):abs(u_y), ...
     [ zeros(length(dof_free),1), stress_hist(2,:) *1e-3 ], 'rx--')
xlabel('Displacement, [mm]')
ylabel('Kirchhoff stress, [kN/mm2]')
