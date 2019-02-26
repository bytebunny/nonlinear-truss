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
[el_x, el_y] = coordxtr(e_dof, coord, dof, 2);
dof_free = [3 4]; % Free DOFs.
bc = [1 0;
      2 0;
      5 0;
      6 0]; % CALFEM: Boundary conditions matrix (DOF, prescribed value).

  % Load and time stepping:
u_y = 2*norm(coord(1,:)-coord(2,:)); % Displace by twice the initial
                                     % length (see fig. 3.8b on p.90).
dt = 0.1; % [s].
total_t = 1; % [s].

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
u_hist = u(dof_free); % Only node 2 can move in X and Y.
bar_forces_history = zeros( 1, size(e_dof,1) ); % Element bar forces.

for t=0:dt:total_t
    du(dof_free) = [0; dt * u_y];
    u_el = extract(e_dof,du); % Element displacements.
    du_el = extract(e_dof,du); % Incremental element displacements.
    
    for iter=1:max_iter % Newton's iteration:
        K_system = sparse(zeros( numel(dof) )); % Global stiffness
                                                % matrix.
        residual_vector = zeros(numel(dof),1);
        for i=1:size(e_dof,1) % Loop over elements:
            [force,K] = element_routine(el_x(i,:), el_y(i,:), ...
                                        u_el(i,:), du_el(i,:), ...
                                        params,i);
keyboard
        end
        
        if error <= error_tol
            fprintf('/// Step %d converged in %d iterations.',...
                    t*total_t/dt,iter)
            break
        end
    end % Newton's iteration.
    
    if iter == max_iter
        error('/// Step %d did NOT converge in %d iterations.',...
              t*total_t/dt,max_iter)
    end
end

%% Post-processing