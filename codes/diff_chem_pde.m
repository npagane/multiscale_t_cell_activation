function [c,f,s] = diff_chem_pde(r,t,u,DuDr,...
    D_IL2, k_endo, K_IL2_tr, k_deg, k_IL2_inj, D_tr_interp, chi_0, K_cc, k_prolif_max, K_prol_IL2, U_min, n_U, r_half,...
    K_JAK_pSTAT5, n_JAK_pSTAT5, K_IL2_JAK_high, K_IL2_JAK_low, R_IL2R_alpha_0,...
     k_transc_basal_mIL2R_alpha_tr, k_transc_pSTAT5_mIL2R_alpha, k_deg_mIL2R_alpha, k_transl_IL2R_alpha, k_deg_IL2R_alpha) %parameters


S_JAK = (u(4)/u(2))/(u(4)/u(2)+R_IL2R_alpha_0)*u(1)/(u(1)+K_IL2_JAK_high*(u(4)/u(2))/(u(4)/u(2)+R_IL2R_alpha_0)+K_IL2_JAK_low*R_IL2R_alpha_0/(u(4)/u(2)+R_IL2R_alpha_0));
P_on_JAK_pSTAT5 = 1/(1+(K_JAK_pSTAT5/S_JAK)^n_JAK_pSTAT5);
c = [1
    1
    u(2)/u(3)
    u(2)/u(4)];

f = [D_IL2*DuDr(1)
    -u(2)*chi_0*(1-u(2)/K_cc)*DuDr(1)+D_tr_interp(r)*DuDr(2)+u(2)*dU_dr(r,n_U,U_min,r_half)*(1-u(2)/K_cc)
    -u(2)*chi_0*(1-u(2)/K_cc)*DuDr(1)+D_tr_interp(r)*DuDr(2)+u(2)*dU_dr(r,n_U,U_min,r_half)*(1-u(2)/K_cc)
    -u(2)*chi_0*(1-u(2)/K_cc)*DuDr(1)+D_tr_interp(r)*DuDr(2)+u(2)*dU_dr(r,n_U,U_min,r_half)*(1-u(2)/K_cc)];

    C_IL2R = 1/2*(u(4)/u(2) + 5000 + 2700) -1/2*sqrt((u(4)/u(2) + 5000 + 2700)^2-4*u(4)/u(2)*5000);  
    s = [-C_IL2R*u(2)*k_endo*u(1)/(K_IL2_tr+u(1))-k_deg*u(1) + k_IL2_inj
     k_prolif_max*P_on_JAK_pSTAT5*u(2)/(K_prol_IL2+P_on_JAK_pSTAT5)*(1-u(2)/K_cc)
    u(2)^2/u(3)*(k_transc_basal_mIL2R_alpha_tr + P_on_JAK_pSTAT5*k_transc_pSTAT5_mIL2R_alpha - k_deg_mIL2R_alpha*(u(3)/u(2))) +  k_prolif_max*P_on_JAK_pSTAT5*u(2)/(K_prol_IL2+P_on_JAK_pSTAT5)*(1-u(2)/K_cc)
    u(2)^2/u(4)*(k_transl_IL2R_alpha*u(3)/u(2) - k_deg_IL2R_alpha*(u(4)/u(2))) +  k_prolif_max*P_on_JAK_pSTAT5*u(2)/(K_prol_IL2+P_on_JAK_pSTAT5)*(1-u(2)/K_cc)
    ];

end


