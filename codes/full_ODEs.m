function dx_dt = full_ODEs(t, x,...
        C_TCR_N, R_CD28,...
        k_on_costim, k_off_costim, k_endo_CTLA4 ,...
        k_transc_basal_mIL2R_alpha, k_transc_TCR_costim_mIL2R_alpha,  k_transc_pSTAT5_mIL2R_alpha, k_deg_mIL2R_alpha, ...
        k_transl_IL2R_alpha, k_deg_IL2R_alpha, R_IL2R_bg, K_IL2_a_bg, k_on_IL2, k_off_IL2, k_endo_C_IL2R, k_transc_TCR_mIL2, ...
        k_deg_mIL2, ... 
        K_TCR_IL2R_alpha, n_TCR_IL2R_alpha, K_JAK_pSTAT5, n_JAK_pSTAT5, K_IL2_JAK_high, K_IL2_JAK_low, R_IL2R_alpha_0, ...
        K_costim_IL2R_alpha, n_costim_IL2R_alpha ,K_TCR_IL2, n_TCR_IL2, K_pSTAT5_IL2, n_pSTAT5_IL2, K_costim_IL2, n_costim_IL2, I,...
        n_tr, n_tr0_trans, V_neighbor, S_DC_TR, f_contact_low, f_contact_high, K_C_CD28_tr_U_min, time_end_priming, C_CD28_tr_interp, C_CTLA4_tr_interp,...
        IL2_pulse_interp, K_TCR_PI3K, n_TCR_PI3K, K_costim_PI3K, n_costim_PI3K, K_JAK_PI3K, n_JAK_PI3K,k_deg_mMyc,k_deg_Myc,k_transc_basal_mMyc,k_transc_PI3K_mMyc,k_transl_Myc)
    
%%%%%%%%%%%%%%%
%%%Variables%%%
%%%%%%%%%%%%%%%
    
%%Costimulation;
S_T = 4*pi*5^2 ;%surface area of T cells, assuming the radious of T cell is 5 um
S_DC = 2000; % surface area of Dendritic cells
L_CD80_86 = x(1);
C_costim =  1/2*(R_CD28/S_T + L_CD80_86/S_DC + k_off_costim/k_on_costim - sqrt((R_CD28/S_T + L_CD80_86/S_DC + k_off_costim/k_on_costim)^2 -4*R_CD28/S_T*L_CD80_86/S_DC));
C_CD28_tr = C_CD28_tr_interp(L_CD80_86/S_DC);
C_CTLA4_tr = C_CTLA4_tr_interp(L_CD80_86/S_DC);   

%%Dynamics of IL2 receptor
mR_IL2R_alpha =  x(2);
R_IL2R_alpha = x(3);
C_IL2R = x(4);

%%Dynamics of IL2 secretion
mI_IL2 = x(5);

S_JAK = R_IL2R_alpha/(R_IL2R_alpha+R_IL2R_alpha_0)*(I+IL2_pulse_interp(t))/(I+IL2_pulse_interp(t)+K_IL2_JAK_high*R_IL2R_alpha/(R_IL2R_alpha+R_IL2R_alpha_0)+K_IL2_JAK_low*R_IL2R_alpha_0/(R_IL2R_alpha+R_IL2R_alpha_0));

%%Dynamics of Myc
mR_Myc = x(6);
R_Myc = x(7);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%activation functions%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%Dynamics of IL2 receptor
P_on_TCR_IL2R_alpha = 1/(1+(K_TCR_IL2R_alpha/C_TCR_N)^n_TCR_IL2R_alpha);
P_on_TCR_costim_IL2R_alpha = 1/(1+((C_TCR_N/K_TCR_IL2R_alpha)^n_TCR_IL2R_alpha+(C_costim/K_costim_IL2R_alpha)^n_costim_IL2R_alpha)^(-1));
if t <= time_end_priming  
    P_on_JAK_pSTAT5 = 1/(1+(K_JAK_pSTAT5/S_JAK)^n_JAK_pSTAT5)*(1-P_on_TCR_IL2R_alpha);
else
    P_on_JAK_pSTAT5 = 1/(1+(K_JAK_pSTAT5/S_JAK)^n_JAK_pSTAT5);
end

%IL2 forced to priming T cells
if false
    P_on_JAK_pSTAT5 = 0.8;
end

%P_on_costim_IL2R_alpha = 1; % turn off the effect of costimulation
%%dynamics of IL2 secretion
P_on_pSTAT5_IL2 = 1/(1+(K_pSTAT5_IL2/double(P_on_JAK_pSTAT5))^n_pSTAT5_IL2);
P_on_TCR_IL2 = 1/(1+(K_TCR_IL2/C_TCR_N)^n_TCR_IL2)*(1-P_on_pSTAT5_IL2);
P_on_costim_IL2 = 1/(1+(K_costim_IL2/C_costim)^n_costim_IL2);
%P_on_costim_IL2 = 1; % turn off the effect of costimulation
P_on_signals_PI3K = ( 1 + ( (C_TCR_N/K_TCR_PI3K)^n_TCR_PI3K + (C_costim/K_costim_PI3K)^n_costim_PI3K + (S_JAK/K_JAK_PI3K)^n_JAK_PI3K )^(-1))^(-1);
%P_on_TCR_PI3K = 1/(1+(K_TCR_PI3K/C_TCR_N)^n_TCR_PI3K);
%P_on_costim_PI3K = 1/(1+(K_costim_PI3K/C_costim)^n_costim_PI3K);
P_on_JAK_PI3K = ( 1 + ( (S_JAK/K_JAK_PI3K)^n_JAK_PI3K )^(-1))^(-1);


%%%%%%%%%%
%%%ODEs%%%
%%%%%%%%%%

%%Costimulation: Linking CD28 the efficiency of transendocytosis.
if n_tr-n_tr0_trans > 0   
    dL_CD80_86_dt = -(f_contact_low+(f_contact_high*f_contact_low-f_contact_low)*1/(1+ (C_CD28_tr/K_C_CD28_tr_U_min)^(-1)))*k_endo_CTLA4*C_CTLA4_tr*V_neighbor*(n_tr-n_tr0_trans)*S_DC_TR;

else
    dL_CD80_86_dt =0;
end


if t <= time_end_priming && C_costim >= 5    
    %%Dynamics of IL2 receptor
    dmR_IL2R_alpha_dt = k_transc_basal_mIL2R_alpha + k_transc_TCR_costim_mIL2R_alpha*P_on_TCR_costim_IL2R_alpha + k_transc_pSTAT5_mIL2R_alpha*double(P_on_JAK_pSTAT5) - k_deg_mIL2R_alpha*mR_IL2R_alpha;
    dR_IL2R_alpha_dt = k_transl_IL2R_alpha*mR_IL2R_alpha - k_deg_IL2R_alpha*R_IL2R_alpha;

    %steady state assumption of IL2R
    R_IL2R = 1/2*(R_IL2R_alpha + R_IL2R_bg + K_IL2_a_bg) -1/2*sqrt((R_IL2R_alpha + R_IL2R_bg + K_IL2_a_bg)^2-4*R_IL2R_alpha*R_IL2R_bg);
    d_C_IL2R_dt = k_on_IL2*(R_IL2R-C_IL2R)*(I+IL2_pulse_interp(t)) -k_off_IL2*C_IL2R- k_endo_C_IL2R*C_IL2R;
    
    %Dynamics of IL2 secretion
    d_mI_IL2_dt = k_transc_TCR_mIL2*P_on_TCR_IL2*P_on_costim_IL2 - k_deg_mIL2*mI_IL2;
    
    % Dynamics of Myc
    dmR_Myc_dt =  k_transc_basal_mMyc + k_transc_PI3K_mMyc*P_on_signals_PI3K - k_deg_mMyc*mR_Myc;
    dR_Myc_dt = k_transl_Myc*mR_Myc - k_deg_Myc*R_Myc;
   
else
    dmR_IL2R_alpha_dt = k_transc_pSTAT5_mIL2R_alpha*double(P_on_JAK_pSTAT5) - k_deg_mIL2R_alpha*mR_IL2R_alpha;
    dR_IL2R_alpha_dt = k_transl_IL2R_alpha*mR_IL2R_alpha - k_deg_IL2R_alpha*R_IL2R_alpha;

    %steady state assumption of IL2R
    R_IL2R = 1/2*(R_IL2R_alpha + R_IL2R_bg + K_IL2_a_bg) -1/2*sqrt((R_IL2R_alpha + R_IL2R_bg + K_IL2_a_bg)^2-4*R_IL2R_alpha*R_IL2R_bg);
    d_C_IL2R_dt = k_on_IL2*(R_IL2R-C_IL2R)*(I+IL2_pulse_interp(t)) -k_off_IL2*C_IL2R- k_endo_C_IL2R*C_IL2R;

    %%Dynamics of IL2 secretion
    d_mI_IL2_dt =  - k_deg_mIL2*mI_IL2;
    
    % Dynamics of Myc
    dmR_Myc_dt =  k_transc_basal_mMyc + k_transc_PI3K_mMyc*P_on_JAK_PI3K - k_deg_mMyc*mR_Myc;
    dR_Myc_dt = k_transl_Myc*mR_Myc - k_deg_Myc*R_Myc;
    
end

dx_dt = double([dL_CD80_86_dt; dmR_IL2R_alpha_dt ; dR_IL2R_alpha_dt ; d_C_IL2R_dt ; d_mI_IL2_dt; dmR_Myc_dt ; dR_Myc_dt]);
        
end

