function [f,k] = element_routine(x, y, u, du, params, el_num)
%% ELEMENT_ROUTINE Returns element bar force and stiffness in 2D.
% Input:
% x - initial X coordinates of nodes.
% y - same but Y.
% u - element displacements.
% du - increment of element displacements.
% params - structure with material and cross-sectional parameters.
% el_num - element number.
%
% Output:
% f - element force vector.
% k - element stiffness matrix.

% Update node positions:
x_new = x + u(1:2:end) + du(1:2:end);
y_new = y + u(2:2:end) + du(2:2:end);
nodes = [x; y];
nodes_new = [x_new; y_new];

diff = nodes(:,2) - nodes(:,1);
diff_new = nodes_new(:,2) - nodes_new(:,1);
l = norm(diff); % Initial length.
l_new = norm(diff_new); % Current length.

n_new = diff_new / l_new; % Element normal vector.

% From Example 3.1 on p.67 in Bonet & Wood:
area_new = params.areas(el_num) * (l_new / l)^(-2*params.nu);

vol = params.areas(el_num) * l;
vol_new = vol * (l_new / l)^(1 - 2 * params.nu); % Eq. (3.4b) in Bonet & Wood.
J = vol_new / vol;

strain = log(l_new / l);
% Compute stress and material stiffness:
[tau, Dtau_Deps] = material_routine(strain, params);

Df_Du = ( area_new / l_new * Dtau_Deps / J - ...
          2 * area_new / l_new * tau / J ) * (n_new * n_new') + ...
        tau / J * area_new / l_new * eye(2);

T = tau * vol / l_new * n_new; % Bar force at node 2.

f = [-T; T];
k = [Df_Du  -Df_Du;
     -Df_Du  Df_Du];
end
