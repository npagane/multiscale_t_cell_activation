function u0 = chempde_D_ic(r, I_tot_interp,   n_tr_tot_interp,n_mIL2Ra_tr_interp, n_IL2Ra_tr_interp)


u0 = [I_tot_interp(r)
     n_tr_tot_interp(r)
     n_mIL2Ra_tr_interp(r)
     n_IL2Ra_tr_interp(r)
     ];

end

