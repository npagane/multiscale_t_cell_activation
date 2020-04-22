%%%TCR signaling via kinetic proofreading.
function F = TCR_func(S_SHP1,L_antigen, R_TCR, kappa, nu, N, phi_max, b, gamma, S_SHP1_t, C_TCR_1_star)
r_p_s = (phi_max + b + gamma*S_SHP1 + nu + sqrt((phi_max + b + gamma*S_SHP1 + nu)^2 - 4*phi_max*(b + gamma*S_SHP1)))/(2*(b+gamma*S_SHP1));
r_n_s = (phi_max + b + gamma*S_SHP1 + nu - sqrt((phi_max + b + gamma*S_SHP1 + nu)^2 - 4*phi_max*(b + gamma*S_SHP1)))/(2*(b+gamma*S_SHP1));
a_n_s = (1-r_n_s)*kappa*R_TCR*L_antigen/(nu+kappa*R_TCR)*(1-(r_n_s/r_p_s)^(N+1))^(-1);
a_p_s = -a_n_s*(r_n_s/r_p_s)^(N+1)*(r_p_s-1)/(r_n_s-1);
C_TCR_1_s = a_p_s*r_p_s + a_n_s*r_n_s;
F = S_SHP1_t*C_TCR_1_s/(C_TCR_1_s+C_TCR_1_star) - S_SHP1;



