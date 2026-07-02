%solve SLHR multistrain
function [eqsolns] = solve_multi_EE(alpha_1, alpha_2, beta_1, beta_2, gamma_1, gamma_2, iota)
syms S L_1 L_2 H_1 H_2;
N_0 = 300;
D = alpha_1 + alpha_2;
R_var_sub = N_0 - S - L_1 - L_2 - H_1 - H_2;

dSdt   = - beta_1*S*H_1 - beta_2*S*H_2 + iota*R_var_sub;
dL_1dt = (beta_1*S*H_1 - alpha_1*L_1)*D - alpha_2*beta_2*L_1*H_2;
dL_2dt = (beta_2*S*H_2 - alpha_2*L_2)*D - alpha_1*beta_1*L_2*H_1;
dH_1dt = (alpha_1*L_1 - gamma_1*H_1)*D + alpha_1*beta_1*L_2*H_1;
dH_2dt = (alpha_2*L_2 - gamma_2*H_2)*D + alpha_2*beta_2*L_1*H_2;

equilibria = solve(dSdt, dL_1dt, dH_2dt, dH_1dt, dL_2dt, S, L_1, L_2, H_1, H_2);
%print(equilibria)
S_eq = equilibria.S;
L_1_eq = equilibria.L_1;
L_2_eq = equilibria.L_2;
%L_12_eq = equilibria.L_12;
H_1_eq = equilibria.H_1;
H_2_eq = equilibria.H_2;
%R_eq = equilibria.R;
eqsolns = [S_eq, L_1_eq, L_2_eq, H_1_eq, H_2_eq];
end