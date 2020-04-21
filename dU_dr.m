%%%Attractive force function
function f = dU_dr(r,n,U_min, r_half)
f = n*U_min*r.^(n-1)./(r_half^n+r.^n)-n*U_min*r.^(2*n-1)./(r_half^n+r.^n).^2;
end
