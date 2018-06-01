function [Sh_cell] = t_sweep(pre_calc_Shock, i0, N_steps, t0, output_file)
% A function that find shocks when sweeping in t. 
% The sweeps are done with quadratic step sizes in t since the only t
% dependence is through sqrt(nustar*t).

%This is to save all the intermediary results to file.
if nargin==4
    output_file='data/t_sweep-tmp.tsv';
end
delete(output_file)

%initializing the output cell
Sh_cell=cell(N_steps+1,1);
Sh_cell{1}=pre_calc_Shock;

%saving initial data
data=[Sh_cell{1}.t, Sh_cell{1}.psimax, Sh_cell{1}.psimin];
save(output_file, 'data', '-ascii','-append')

%Retrieveing the constant shock properties
Z=pre_calc_Shock.Z; n=pre_calc_Shock.n; m=pre_calc_Shock.m;
Mach=pre_calc_Shock.M; taui=pre_calc_Shock.taui;
nu_star=pre_calc_Shock.nu_star; tol=pre_calc_Shock.tol;


for j=1:N_steps
    psiminmax_prev=[Sh_cell{j}.psimax, Sh_cell{j}.psimin];%New initial guess
    t=t0*(j+i0)^2; %step in t
    %Printout to show progress
    fprintf('============================================================\n')
    fprintf('j+i0 = %d\t t = %1.2f\n\n',j+i0, t)
    %Calc new shock
    Sh_cell{j+1}=Shock_pkg_new.Shock_col(Z,n,m, taui,Mach, t, nu_star, psiminmax_prev, tol);
    %saves the new data to file
    data=[t, Sh_cell{j+1}.psimax, Sh_cell{j+1}.psimin];
    save(output_file, 'data', '-ascii','-append')
end
end

