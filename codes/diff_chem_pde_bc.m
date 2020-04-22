function [ pl,ql,pr,qr ] = chempde_D_bc(rl,ul,rr,ur,tspan_diff, R_IL2R_tc, C_IL2R_tc, K_IL2 ,q_IL2, D_IL2 ,n_tr0)




pl = [q_IL2 - 1.83*10^2/3600*ul(1)*(R_IL2R_tc-C_IL2R_tc) + 0.83/3600*C_IL2R_tc% - 1.83*10^2/3600*ul(1)*(R_IL2R_tc-C_IL2R_tc) + 0.83/3600*C_IL2R_tc %IL2
      0 %Tregs
      0 %mIL2Ra
      0 %IL2Ra
      ]; 
ql = [4*pi*5^2 %IL2
     1 %Tregs
     1 %mIL2Ra
     1 %IL2Ra
     ];

pr = [0; 0; 0; 0];
qr = [1; 1; 1; 1];


end

