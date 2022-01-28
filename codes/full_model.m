%%% These are Matlab codes for the model of T cell activations via three
%%% cell- interactions. Among various internal version, this is the version
%%% used in the publication "A local regulatory T cell feedback circuit 
%%% maintains immunological homeostasis by pruning self-activated T cells". 
%%% For detailed mathematical descriptions and biological assumptions, 
%%% please refer to the the paper.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Loading parameter combinations%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('prmset_off_corrected.mat'); %load partameter configurations
idx_prm = 938; % choose a parameter configuration (1~20000)
disp(strcat(prmset{idx_prm,1}));

%%%%%%%%%%%%%%%%
%%%parameters%%%
%%%%%%%%%%%%%%%%
%%TCR signaling
L_antigen = prmset{idx_prm, 17};
R_TCR = 30000;
kappa = 10^(-4); %2D: 1.46*10^(-2);
nu = prmset{idx_prm, 2};
N = 5;
phi_max = 0.09;
b = 0.04;
gamma = 1.2*10^(-6);
S_SHP1_t = 60000;
C_TCR_1_star = 500;

%%Costimulation signaling
R_CD28 = 30000;
L_CD80_86_i = scaled_inflam*prmset{idx_prm, 3}; %CD86
R_CTLA4_max = 24500;
k_on_costim = 0.77; %1/s
k_off_costim = 28;%1/s
k_on_CTLA4 = 1.09;%1/s
k_off_CTLA4 = 5.1; %0.0082; %1/s
k_syn_CTLA4 = 12830; %1/m
k_endo_CTLA4 = 0.291*60;%1/h
k_deg_CTLA4 = 0.0049;%1/m
k_recyc_CTLA4 = 0.021;%1/m
V_neighbor = 4/3*pi*(15^3-5^3);%um^3
S_DC_TR = 8; %average contact area between DCs and Tregs. um^2
f_contact_low = prmset{idx_prm, 4};
f_contact_high = prmset{idx_prm, 5};

%%Dynamics of IL2 receptor
k_transc_basal_mIL2R_alpha = 0.0015;
k_transc_basal_mIL2R_alpha_tr = 0.5;
k_transc_pSTAT5_mIL2R_alpha = 20; 
k_transc_TCR_costim_mIL2R_alpha = 0.5; % OR gate of TCR and costimulation
k_deg_mIL2R_alpha = 0.2;
k_transl_IL2R_alpha = 20;
k_deg_IL2R_alpha = 0.05;
R_IL2R_bg = 5000; % assuming constant value
K_IL2_a_bg = 2700;
k_on_IL2 = 183;
k_off_IL2 = 0.83;
k_endo_C_IL2R = 2;

%%dynamics of IL2 secretion
k_transc_TCR_mIL2 = 162*(1+prmset{idx_prm, 11}/66.518);  %compensate
k_transc_costim_mIL2 = 0; % not used since TCR and costimulation are merged as 'AND' gate.
k_deg_mIL2 = 0.2;
k_transl_IL2 = 266;

%%dynamics of PI3K and MYC
k_deg_mMyc = 1/5;
k_deg_Myc = 0.099;
k_transc_basal_mMyc = 0.03*k_deg_Myc;
k_transc_PI3K_mMyc = 50*k_deg_Myc;
k_transl_Myc = 20*k_deg_mMyc;
K_TCR_PI3K = prmset{idx_prm, 12}; %NP EDIT deduced the following params (0.3)
K_costim_PI3K = prmset{idx_prm, 14};
K_JAK_PI3K = prmset{idx_prm, 13}; 
n_TCR_PI3K = 1;
n_costim_PI3K = 1;
n_JAK_PI3K = 1;

%%dynamics of CD86 (not used in this version)
k_transc_mCD80_86 = 20;
k_deg_mCD80_86 = 1/5;
k_transl_CD80_86 = 150*2;
k_deg_CD80_86 = 0.05*2; %to make the initial state as a stationary state


%%%%%%%%%%%%%%%%%%%%%%%%
%%%Upstream signaling%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%Dynamics of IL2 receptor
K_TCR_IL2R_alpha = prmset{idx_prm, 6};
n_TCR_IL2R_alpha = 1;

K_JAK_pSTAT5 = prmset{idx_prm, 7};
n_JAK_pSTAT5 = 2;

K_IL2_JAK_high = 6.02*10^(-4); % molecules per um^3 = 1 pM
K_IL2_JAK_low = 50*6.02*10^(-4); % molecules per um^3 = 50 pM
R_IL2R_alpha_0 = 2*10^4;
K_costim_IL2R_alpha = prmset{idx_prm, 8}; % not yet determined
n_costim_IL2R_alpha = 1;

%%dynamics of IL2 secretion
K_TCR_IL2 = prmset{idx_prm, 9};
n_TCR_IL2 = 1;
K_pSTAT5_IL2 = prmset{idx_prm, 10}; 
n_pSTAT5_IL2 = 1; 
K_costim_IL2 = prmset{idx_prm, 11}; 
n_costim_IL2 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%IL2 diffusion & Tregs chemotaxis%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 2; % spherical coordinates
r_max = 100;
r = linspace(5, r_max,40);

D_IL2 = prmset{idx_prm, 16};
n_tr0_low = prmset{idx_prm, 15}*0.6;
n_tr0_high = 10*prmset{idx_prm, 15}*0.6;
n_tr0 = [n_tr0_high n_tr0_high n_tr0_high repmat(n_tr0_low,1,37)];
n_tr0_trans = 0; % threshold of treg density above which effective transendocytosis happens

k_endo_C_IL2R_tr = 2/3600; % s^-1
R_tr = 10000;
K_IL2_tr = 0.006;% molecules/um^3 
k_deg = 0.1; %h^-1

k_IL2_inj = 0;
K_IL2_tc = 0.006;

%Tregs chemotaxis
D_tr_low = 0.17/10;
D_tr_high = 0.17;
D_tr = [repmat(D_tr_low,1,3) repmat(D_tr_high,1,37)];
D_tr_interp = griddedInterpolant(r,D_tr,'pchip');

chi_0 = 0;
K_cc = 0.003; %carrying capacity ~0.003 /um^3  

%Treg proliferation
k_prolif_max = log(2)/18/3600; %/s, doubling time 18 hours.
K_prol_IL2 = 0.1; % pSTAT5 level

%Potential near DC
U_min = 0.05;
K_C_CD28_tr_U_min = 5;
K_pSTAT5_tr_U_min = 0.1;
n_U = 5;
r_half = 15;

%%Treg internal states
%initial condition 
R_mIL2Ra_tr0 = 2.5; % steady state
R_IL2Ra_tr0 = 1000; % steady state

%%Costimulation;
S_T = 4*pi*5^2 ;%surface area of T cells, assuming the radious of T cell is 5 um
S_DC = 2000; % surface area of Dendritic cells


%%%%%%%%%%%%%%%%%%%%%%
%% Simulation setup%%%
%%%%%%%%%%%%%%%%%%%%%%

time_i = 0; %in hours
time_f = 120;%in hours
time_end_priming = 48; %in hours
time_del = 5/60; % time interval for seperation of time scales
tspan_diff = linspace(0,time_del*3600,3); %for chemo_diffusion equation
tspan_tot = [time_i time_f];
num_t_step = (time_f-time_i)/time_del;

t_tot_diff = linspace(time_i,time_f,num_t_step+1);
n_tr_tot = zeros(num_t_step+1,40);
n_tr_tot(1,:) = n_tr0;
D_tr_low_tot = zeros(num_t_step+1,1);
D_tr_low_tot(1) = D_tr_low;
n_mIL2Ra_tr_tot = zeros(num_t_step+1,40);
n_mIL2Ra_tr_tot(1,:) = R_mIL2Ra_tr0*n_tr_tot(1,:);
n_IL2Ra_tr_tot = zeros(num_t_step+1,40);
n_IL2Ra_tr_tot(1,:) = R_IL2Ra_tr0*n_tr_tot(1,:);

x_tot = zeros(200000,7);  
t_tot = zeros(200000,1);
I_tot = zeros(num_t_step+1,40);
I_tot(1,:) = 0;

x_tot(1,:) = [L_CD80_86_i 0 0 0 0 0 0]; % initial condition
t_tot(1) = 0;

%count the number of time points
num_pt = 1;

%%dynamics of IL2 pulse NP EDIT HERE
%velocity = 20.0; 
vlow = 5; vmax = 20; coefficient = 0.6;
velocity = coefficient*vlow + (1-coefficient)*vmax;
freq_IL2_pulse = 6.2 *(velocity/vlow); %per hour
duration_IL2_pulse = 2 / 60; %(10/((1-coefficient)*vlow)) / 60; % minutes (but in hours)

time_pulse = time_end_priming;
trajectory_pulse = [time_pulse];
signal_pulse = [0];
IL2_signal_pulse = 10 * 0.000602;
eps = 1e-10;
shape = 1.6750;
scale = 2.0786;
e50il2 = 30 *  0.000602;
hill_n = 1;

% while time_pulse < time_f
%     time_pulse = time_pulse + 1/freq_IL2_pulse;
%     trajectory_pulse = [trajectory_pulse, time_pulse-eps]; signal_pulse = [signal_pulse,0];
%     trajectory_pulse = [trajectory_pulse, time_pulse+eps]; signal_pulse = [signal_pulse,IL2_signal_pulse];
%     time_pulse = time_pulse + duration_IL2_pulse;
%     trajectory_pulse = [trajectory_pulse, time_pulse-eps]; signal_pulse = [signal_pulse,IL2_signal_pulse];
%     trajectory_pulse = [trajectory_pulse, time_pulse+eps]; signal_pulse = [signal_pulse,0];
% end

ld = 10;
signal_temp = 0;
integrated_signal = 0;
while time_pulse < time_f
    coefficient = 1-(integrated_signal(end)^hill_n)/(e50il2^hill_n+integrated_signal(end)^hill_n);
    velocity = coefficient*vlow + (1-coefficient)*vmax;
    memory_time = ld/velocity; % memory time inversely dependenton frequency
    freq_IL2_pulse = 6.2 *(velocity/vlow); 
    t_temp = exprnd(freq_IL2_pulse);
    t_temp = 1/t_temp;
    duration_IL2_pulse = gamrnd(shape,scale,1)/60;
    if t_temp > duration_IL2_pulse
        trajectory_pulse = [trajectory_pulse, time_pulse+duration_IL2_pulse-eps]; signal_pulse = [signal_pulse,signal_temp];
        integrated_signal = integrated_signal+signal_temp;
        signal_temp = 0;
        trajectory_pulse = [trajectory_pulse, time_pulse+duration_IL2_pulse+eps]; signal_pulse = [signal_pulse,signal_temp];
        integrated_signal = max(0,integrated_signal - IL2_signal_pulse*(t_temp/memory_time));
        time_pulse = time_pulse + t_temp;
        trajectory_pulse = [trajectory_pulse, time_pulse-eps]; signal_pulse = [signal_pulse,signal_temp];
        integrated_signal = integrated_signal+signal_temp;
        signal_temp = IL2_signal_pulse;
        trajectory_pulse = [trajectory_pulse, time_pulse+eps]; signal_pulse = [signal_pulse,signal_temp];
        integrated_signal = integrated_signal+signal_temp;
    else
        time_pulse = time_pulse + t_temp;
        trajectory_pulse = [trajectory_pulse, time_pulse-eps]; signal_pulse = [signal_pulse,signal_temp];
        integrated_signal = integrated_signal+signal_temp;
        signal_temp = signal_temp + IL2_signal_pulse;
        trajectory_pulse = [trajectory_pulse, time_pulse+eps]; signal_pulse = [signal_pulse,signal_temp];
        integrated_signal = integrated_signal+signal_temp;
    end
end

%plot(trajectory_pulse, signal_pulse)
IL2_pulse_interp = griddedInterpolant(trajectory_pulse, signal_pulse, 'pchip');
%integrated_interp = griddedInterpolant(trajectory_pulse, integrated_signal, 'pchip');
%plot(trajectory_pulse, IL2_pulse_interp(trajectory_pulse))
%axis([0 time_f 0 2e-2])
%figure; plot(trajectory_pulse, integrated_interp(trajectory_pulse))
%IL2_pulse_interp = griddedInterpolant([0,time_f],[0,0],'pchip');

%%%%%%%%%%%%%%%%
%% Simulation%%%
%%%%%%%%%%%%%%%%

%%Solve for TCR signaling 
fun = @(x) TCR_func(x, L_antigen, R_TCR, kappa, nu, N, phi_max, b, gamma, S_SHP1_t, C_TCR_1_star);
x_seed = 400000;
options = optimset('Display','off');
S_SHP1 = fsolve(fun, x_seed,options);
r_p = (phi_max + b + gamma*S_SHP1 + nu + sqrt((phi_max + b + gamma*S_SHP1 + nu)^2 - 4*phi_max*(b + gamma*S_SHP1)))/(2*(b+gamma*S_SHP1));
r_n = (phi_max + b + gamma*S_SHP1 + nu - sqrt((phi_max + b + gamma*S_SHP1 + nu)^2 - 4*phi_max*(b + gamma*S_SHP1)))/(2*(b+gamma*S_SHP1));
a_n = (1-r_n)*kappa*R_TCR*L_antigen/(nu+kappa*R_TCR)*(1-(r_n/r_p)^(N+1))^(-1);
a_p = -a_n*(r_n/r_p)^(N+1)*(r_p-1)/(r_n-1);
C_TCR_N = a_p*r_p^N + a_n*r_n^N;


%%Build interploted functions for C_CTLA4_tr and C_CD28_tr based on initial
%%L_CD80_86 to reduce the computation time.
L_CD80_86_surf = 0:1:(L_CD80_86_i/S_DC+1);
R_CD28_surf_tr =  R_CD28/S_T;
R_CTLA4_surf_tr = R_CTLA4_max/S_T;
costim_complex_sol = zeros(length(L_CD80_86_surf),3);
costim_complex_sol(:,1) = L_CD80_86_surf;
for L_CD80_86_surf_idx = L_CD80_86_surf
  fun = @(x) CD28_CTLA4_func(x,  R_CD28_surf_tr, R_CTLA4_surf_tr, L_CD80_86_surf_idx,...
    k_on_costim, k_off_costim, k_on_CTLA4, k_off_CTLA4);
    x_seed = [10,10];
    %temp = fsolve(fun, x0,options);
    temp = fsolve(fun, x_seed, options);
    costim_complex_sol(L_CD80_86_surf_idx+1,2:3)= temp;
end
C_CD28_tr_interp = griddedInterpolant(costim_complex_sol(:,1),costim_complex_sol(:,2),'pchip');
C_CTLA4_tr_interp = griddedInterpolant(costim_complex_sol(:,1),costim_complex_sol(:,3),'pchip');

%%ODE and PDE solvers
for i = 1:num_t_step
    tspan = [(i-1)*time_del i*time_del]; %in hours
  
    %Initial conditions
    x0 = x_tot(num_pt,:); 
    I_0 = I_tot(i, 1);
    n_tr_l = (n_tr_tot(i,1) + n_tr_tot(i,2) + n_tr_tot(i,3))/3; %average in 5~15um
   
    n_tr_tot_interp = griddedInterpolant(double(r),double(n_tr_tot(i,:)),'nearest');
    I_tot_interp = griddedInterpolant(double(r),double(I_tot(i,:)),'nearest');
    n_mIL2Ra_tr_interp = griddedInterpolant(double(r),double(n_mIL2Ra_tr_tot(i,:)),'nearest');
    n_IL2Ra_tr_interp = griddedInterpolant(double(r),double(n_IL2Ra_tr_tot(i,:)),'nearest');
    
    
    R_IL2Ra_tr_temp = n_IL2Ra_tr_interp(r)/n_tr_tot_interp(r);
    S_JAK_tr_temp = R_IL2Ra_tr_temp./(R_IL2Ra_tr_temp+R_IL2R_alpha_0).*I_tot_interp(r)./(I_tot_interp(r)+K_IL2_JAK_high.*(R_IL2Ra_tr_temp)./(R_IL2Ra_tr_temp+R_IL2R_alpha_0)+K_IL2_JAK_low*R_IL2R_alpha_0./(R_IL2Ra_tr_temp+R_IL2R_alpha_0));
    P_on_TCR_IL2R_alpha_tr_temp = 0;
    P_on_JAK_pSTAT5_tr_temp = 1./(1+(K_JAK_pSTAT5./S_JAK_tr_temp).^n_JAK_pSTAT5)*(1-P_on_TCR_IL2R_alpha_tr_temp);
    P_on_JAK_pSTAT5_tr_l = (P_on_JAK_pSTAT5_tr_temp(1) +P_on_JAK_pSTAT5_tr_temp(2)+P_on_JAK_pSTAT5_tr_temp(3))/3; % average in 5~15um
    
    %%ODEs
    [t,x] = ode45( @(tspan, x)full_ODEs(tspan, x,...
        C_TCR_N, R_CD28,...
        k_on_costim, k_off_costim, k_endo_CTLA4, ...
        k_transc_basal_mIL2R_alpha, k_transc_TCR_costim_mIL2R_alpha,  k_transc_pSTAT5_mIL2R_alpha, k_deg_mIL2R_alpha, ...
        k_transl_IL2R_alpha, k_deg_IL2R_alpha, R_IL2R_bg, K_IL2_a_bg, k_on_IL2, k_off_IL2, k_endo_C_IL2R, k_transc_TCR_mIL2, ...
        k_deg_mIL2, ... 
        K_TCR_IL2R_alpha, n_TCR_IL2R_alpha, K_JAK_pSTAT5, n_JAK_pSTAT5, K_IL2_JAK_high, K_IL2_JAK_low, R_IL2R_alpha_0, ...
        K_costim_IL2R_alpha, n_costim_IL2R_alpha ,K_TCR_IL2, n_TCR_IL2, K_pSTAT5_IL2, n_pSTAT5_IL2, K_costim_IL2, n_costim_IL2, I_0,...
        n_tr_l, n_tr0_trans, V_neighbor, S_DC_TR, f_contact_low, f_contact_high, K_C_CD28_tr_U_min, time_end_priming , C_CD28_tr_interp, C_CTLA4_tr_interp,...
        IL2_pulse_interp, K_TCR_PI3K, n_TCR_PI3K, K_costim_PI3K, n_costim_PI3K, K_JAK_PI3K, n_JAK_PI3K,k_deg_mMyc,k_deg_Myc,k_transc_basal_mMyc,k_transc_PI3K_mMyc,k_transl_Myc),...
        tspan, x0);
     
    %%PDEs    
    
    q_IL2 = k_transl_IL2*x(end,5)/3600;
    R_IL2R_tc = 1/2*(x(end,3) + R_IL2R_bg + K_IL2_a_bg) -1/2*sqrt((x(end,3) + R_IL2R_bg + K_IL2_a_bg)^2-4*x(end,3)*R_IL2R_bg);
    C_IL2R_tc = x(end,4);
    C_CD28_tr = C_CD28_tr_interp(x(end,1)/S_DC);
    
    D_tr_low_tot(i) = D_tr_high -(D_tr_high-D_tr_low)*1/(1+ (C_CD28_tr/K_C_CD28_tr_U_min + P_on_JAK_pSTAT5_tr_l/K_pSTAT5_tr_U_min)^(-1));
    D_tr = [repmat(D_tr_low_tot(i),1,3) repmat(D_tr_high,1,37)];
    D_tr_interp = griddedInterpolant(r,D_tr,'pchip');
    
    sol = pdepe(m,@(r,tspan_diff,u,DuDr) diff_chem_pde( r,tspan_diff,u,DuDr, D_IL2, k_endo_C_IL2R_tr, K_IL2_tr, k_deg, k_IL2_inj, D_tr_interp, chi_0, K_cc,k_prolif_max, K_prol_IL2, U_min*1/(1+ (C_CD28_tr/K_C_CD28_tr_U_min + P_on_JAK_pSTAT5_tr_l/K_pSTAT5_tr_U_min)^(-1)),n_U,r_half,...
        K_JAK_pSTAT5, n_JAK_pSTAT5, K_IL2_JAK_high, K_IL2_JAK_low, R_IL2R_alpha_0,...
        k_transc_basal_mIL2R_alpha_tr/3600, k_transc_pSTAT5_mIL2R_alpha/3600, k_deg_mIL2R_alpha/3600, k_transl_IL2R_alpha/3600, k_deg_IL2R_alpha/3600),... 
        @(r) diff_chem_pde_ic(r, I_tot_interp,  n_tr_tot_interp,n_mIL2Ra_tr_interp, n_IL2Ra_tr_interp ),...
        @(rl,Il,rr,Ir,tspan_diff) diff_chem_pde_bc(rl,Il,rr,Ir,tspan_diff, R_IL2R_tc, C_IL2R_tc ,K_IL2_tc, q_IL2, D_IL2, n_tr0),...
        r, tspan_diff);

    I = sol(:,:,1);
    n_tr = sol(:,:,2);
    n_mIL2Ra_tr = sol(:,:,3);
    n_IL2Ra_tr = sol(:,:,4);
      
    %ODE
    x_tot((num_pt+1):(num_pt+ size(t,1)-1),:) = x(2:end,:);
    t_tot((num_pt+1):(num_pt+ size(t,1)-1)) = t(2:end);
    
    %PDE
    I_tot(i+1,:) = I(3,:);
    n_tr_tot(i+1,:) = n_tr(3,:);
    n_mIL2Ra_tr_tot(i+1,:) = n_mIL2Ra_tr(3,:);
    n_IL2Ra_tr_tot(i+1,:) = n_IL2Ra_tr(3,:);
    
    num_pt = num_pt + size(t,1)-1;
    
    disp(i)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Combine and save simulation output%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_tot = x_tot(1:num_pt,:);
t_tot = t_tot(1:num_pt);
x_tot = [x_tot(:,1:7), k_transl_IL2*x_tot(:,5)/3600]; %IL2 secretion (s^-1)

%pSTAT5 activation of Tc
I_tot_1 = t_tot;
idx = zeros(num_t_step+1,1);
%idx(1) = 1;
for i = 1:num_t_step
    idx(i+1) = find(t_tot == i*time_del);
    I_tot_1((idx(i)+1):idx(i+1))= I_tot(i,1) + IL2_pulse_interp(i*time_del);
    
end

S_JAK = x_tot(:,3)./(x_tot(:,3)+R_IL2R_alpha_0).*I_tot_1./(I_tot_1+K_IL2_JAK_high*x_tot(:,3)./(x_tot(:,3)+R_IL2R_alpha_0)+K_IL2_JAK_low*R_IL2R_alpha_0./(x_tot(:,3)+R_IL2R_alpha_0));
P_on_TCR_IL2R_alpha = 1/(1+(K_TCR_IL2R_alpha/C_TCR_N)^n_TCR_IL2R_alpha);
C_costim = 1/2*(R_CD28/S_T + x_tot(:,1)/S_DC + k_off_costim/k_on_costim - sqrt((R_CD28/S_T + x_tot(:,1)/S_DC + k_off_costim/k_on_costim).^2 -4*R_CD28/S_T*x_tot(:,1)/S_DC));
t_tot_diseng = t_tot(C_costim < 5);
if size(t_tot_diseng,1) ~=0
    if t_tot_diseng(1) < time_end_priming
        time_end_priming = t_tot_diseng(1);
    end
end
P_on_JAK_pSTAT5 = 1./(1+(K_JAK_pSTAT5./S_JAK).^n_JAK_pSTAT5).*(1-P_on_TCR_IL2R_alpha.*(t_tot < time_end_priming));

%%pSTAT5 activation of Tr
R_IL2Ra_tr_temp = n_IL2Ra_tr_tot./n_tr_tot;
S_JAK_tr = R_IL2Ra_tr_temp./(R_IL2Ra_tr_temp+R_IL2R_alpha_0).*I_tot./(I_tot+K_IL2_JAK_high.*(R_IL2Ra_tr_temp)./(R_IL2Ra_tr_temp+R_IL2R_alpha_0)+K_IL2_JAK_low*R_IL2R_alpha_0./(R_IL2Ra_tr_temp+R_IL2R_alpha_0));
P_on_TCR_IL2R_alpha_tr = 0;
P_on_JAK_pSTAT5_tr = 1./(1+(K_JAK_pSTAT5./S_JAK_tr).^n_JAK_pSTAT5)*(1-P_on_TCR_IL2R_alpha_tr);


C_CTLA4 =C_CTLA4_tr_interp(x_tot(:,1)/S_DC);
C_CTLA4_prune = C_CTLA4([1 idx(2:end)']);
C_costim_tr = C_CD28_tr_interp(x_tot(:,1)/S_DC);

IL2_pulse = IL2_pulse_interp(t_tot);
%%Save data
%save('output.mat',...
%        'n_tr_tot', 'D_tr_low_tot', 'n_mIL2Ra_tr_tot', 'n_IL2Ra_tr_tot', 'I_tot', 'x_tot', 't_tot', 't_tot_diff', 'r',...
%        'costim_complex_sol', 'idx', 'I_tot_1', 'S_JAK', 'P_on_TCR_IL2R_alpha', 'C_costim', 'time_end_priming', 'P_on_JAK_pSTAT5', 'S_JAK_tr', 'P_on_JAK_pSTAT5_tr', 'C_CTLA4', 'C_CTLA4_prune', 'C_costim_tr');

%%Plot

figure;
plot(t_tot,x_tot(:,7),'-o');
legend('Myc');
axis([time_i time_f 0 1000]);

if false
    figure;
    subplot(3,1,1);
    plot(t_tot,x_tot(:,1),'-o');
    legend('CD80/CD86');
    %axis([time_i time_f 0 L_CD80_86_i*1.5]);
    subplot(3,1,2);
    plot(t_tot,C_costim,'-o',t_tot, C_costim_tr ,t_tot,C_CTLA4,'-*');
    %axis([time_i time_f 0 100]);
    legend('Costim (TC)','Costim (TR)','CTLA4 (TR)');
    subplot(3,1,3);
    plot(t_tot,x_tot(:,2),'-o',t_tot,x_tot(:,3),'-*',t_tot,x_tot(:,4),'-+');
    legend('mIL2Ra','RIL2Ra', 'CIL2Ra');

    figure;
    plot(t_tot,x_tot(:,6),'-o',t_tot,x_tot(:,7),'-*');
    legend('mMyc','RMyc');
    %axis([time_i time_f 0 1000]);

    figure;
    subplot(3,1,1);
    plot(t_tot,x_tot(:,5),'-.');
    legend('mIL2');
    subplot(3,1,2);
    plot(t_tot,x_tot(:,8),'-s');
    legend('IL2 secretion (molecules/sec)');
    subplot(3,1,3);
    plot(t_tot_diff,I_tot(:,1)/0.000602,'-s');
    legend('IL2 conc (pM)');

    figure;
    subplot(2,1,1);
    plot(t_tot,P_on_JAK_pSTAT5,'-s');
    legend('pSTAT5 of Tc');
    %axis([time_i time_f 0 1]);
    subplot(2,1,2);
    plot(t_tot_diff,P_on_JAK_pSTAT5_tr(:,1),'-o',...
        t_tot_diff,P_on_JAK_pSTAT5_tr(:,2),'-*',t_tot_diff,P_on_JAK_pSTAT5_tr(:,3),'-+');
    %axis([time_i time_f 0 1]);
    legend('pSTAT5 of Tregs at 5um', 'pSTAT5 of Tregs at 10um', 'pSTAT5 of Tregs at 15um');

    figure;
    surf(r,t_tot_diff(1:10:end),I_tot(1:10:end,:)/0.000602);
    ylabel('Time (Hours)'); xlabel('r (um)'); zlabel('Concentration (pM)')

    figure;
    surf(r,t_tot_diff(1:10:end),n_tr_tot(1:10:end,:));
    %axis([0 r_max time_i time_f 0 0.0015])
    ylabel('Time (Hours)'); xlabel('r (um)'); zlabel('Density (#/um^3)')   

    figure
    surf(r,t_tot_diff(1:10:end),n_mIL2Ra_tr_tot(1:10:end,:)./n_tr_tot(1:10:end,:))
    %axis([0 r_max time_i time_f])
    ylabel('Time (Hours)'); xlabel('r (um)'); zlabel('mIL2Ra per cell (tr) (#/cell)')   
    
    figure
    surf(r,t_tot_diff(1:10:end),n_IL2Ra_tr_tot(1:10:end,:)./n_tr_tot(1:10:end,:))
    ylabel('Time (Hours)'); xlabel('r (um)'); zlabel('IL2Ra per cell (tr) (#/cell)') 
    %axis([0 r_max time_i time_f])
    
    figure
    surf(r,t_tot_diff(1:10:end),P_on_JAK_pSTAT5_tr(1:10:end,:))
    ylabel('Time (Hours)'); xlabel('r (um)'); zlabel('pSTAT5 (tr) (a.u.)') 
    %axis([0 r_max time_i time_f])
end
