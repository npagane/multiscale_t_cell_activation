%%%Compe
function F = CD28_CTLA4_func(x, R_CD28_tot_tr, R_CTLA4_tot_tr, L_CD80_86_tot,...
    k_on_costim, k_off_costim, k_on_CTLA4, k_off_CTLA4)
F(1) =  k_on_costim*(R_CD28_tot_tr-x(1))*(L_CD80_86_tot-x(1)-x(2))- k_off_costim*x(1);
F(2) = k_on_CTLA4*(R_CTLA4_tot_tr-x(2))*(L_CD80_86_tot-x(1)-x(2))- k_off_CTLA4*x(2);
