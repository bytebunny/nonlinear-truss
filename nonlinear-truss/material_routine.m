function [stress, mat_stiff] = material_routine(strain, params)
%% MATERIAL_ROUTINE Returns stress and material stiffness.
% Inputs:
% strain
% params - structure containing material parameters.
%
% Output:
% stress
% mat_stiffness - derivative of stress w.r.t. strain.

stress = params.young * strain;
mat_stiff = params.young;
end