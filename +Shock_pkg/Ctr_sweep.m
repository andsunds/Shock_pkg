function [Sh_cell] = Ctr_sweep(pre_calc_shock, N_steps, C_final)

import Shock_pkg.*

Sh_cell=cell(N_steps,1);
Sh_cell{1}=pre_calc_shock;

C0=Sh_cell{1}.trapping_coef;
if nargin==2
    C_final=1;
end
dC=(1-C0)/(N_steps-C_final);

for j=2:N_steps
    Sh_prev=Sh_cell{j-1};
    C_tr=Sh_prev.trapping_coef+dC;
    
    fprintf('============================================================\n')
    fprintf('j = %d,  \t C_tr = %1.3f\n',j,C_tr) %Print out progress
    Sh_cell{j}=Shock_tr(Sh_prev.m,Sh_prev.Z,Sh_prev.n,Sh_prev.tau,...
        C_tr, 'Mach',Sh_prev.Mach, [Sh_prev.phimax, Sh_prev.phimin], Sh_prev.tol);
end

end

