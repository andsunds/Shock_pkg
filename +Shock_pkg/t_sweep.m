function [Sh_cell] = t_sweep(pre_calc_Shock, i0, N_steps, t0, output_file)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin==4
    output_file='data/t_sweep-tmp.tsv';
end
delete(output_file)

Sh_cell=cell(N_steps+1,1);
Sh_cell{1}=pre_calc_Shock;
data=[Sh_cell{1}.t, Sh_cell{1}.phimax, Sh_cell{1}.phimin];
save(output_file, 'data', '-ascii','-append')


m=pre_calc_Shock.m; Z=pre_calc_Shock.Z; n=pre_calc_Shock.n;
Mach=pre_calc_Shock.Mach; tau=pre_calc_Shock.tau;
nu_star=pre_calc_Shock.nu_star; tol=pre_calc_Shock.tol;


for j=1:N_steps
    phiminmax_prev=[Sh_cell{j}.phimax, Sh_cell{j}.phimin];
    t=t0*(j+i0)^2;
    fprintf('============================================================\n')
    fprintf('j+i0 = %d\t t = %1.2f\n\n',j+i0, t)
    Sh_cell{j+1}=Shock_pkg.Shock_col(m,Z,n, tau, 'Mach',Mach, t, nu_star, phiminmax_prev, tol);
    data=[t, Sh_cell{j+1}.phimax, Sh_cell{j+1}.phimin];
    save(output_file, 'data', '-ascii','-append')
end
end

