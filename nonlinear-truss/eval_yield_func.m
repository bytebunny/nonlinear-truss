function [f] = evla_yield_func(tau, tau_y0, H, acc_pl_strain)
%% EVAL_YIELD_FUNC Returns value of the yield function.
% Inputs:
% tau - Kirchhoff stress.
% tau_y0 - initial yield stress.
% H - plastic modulus.
% acc_pl_strain - hardening parameter.
%
% Output:
% f - value of the yield function.

f = abs(tau) - ( tau_y0 + H * acc_pl_strain );
end