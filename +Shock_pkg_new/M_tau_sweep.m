function [Sh_cell] = M_tau_sweep(pre_calc_Shock, iT0,TT, iM0,M0,dM,NM)


Sh_handle=@Shock_pkg_new.Shock_col;


%initializing the cell list to be returned
NT=length(TT);
Sh_cell=cell(NT,NM);
Sh_cell{iT0,iM0}=pre_calc_Shock;

%Retrieveing the constant shock properties
Z=pre_calc_Shock.Z; n=pre_calc_Shock.n; m=pre_calc_Shock.m; t=pre_calc_Shock.t;
nu_star=pre_calc_Shock.nu_star; tol=pre_calc_Shock.tol;

%Printout to show progress
fprintf('============================================================\n')
fprintf('j = %d,  \t tau = %1.3f\n',iT0,TT(iT0)) %Print out progress

Sh_tmp=Shock_pkg_new.Mach_sweep(Sh_handle, Sh_cell{iT0,iM0}, iM0, NM, M0, dM);
for k=1:NM
    Sh_cell{iT0,k}=Sh_tmp{k};
end

for i=(iT0+1):NT
    %Printout to show progress
    fprintf('============================================================\n')
    fprintf('j = %d,  \t tau = %1.3f\n',i,TT(i)) %Print out progress
    try
        psimaxmin_in=[Sh_cell{i-1,iM0}.psimax,Sh_cell{i-1,iM0}.psimin];
        Sh_cell{i,iM0}=Sh_handle(Z,n,m, TT(i), M0, t, nu_star, psimaxmin_in, tol);
        Sh_tmp=Shock_pkg_new.Mach_sweep(Sh_handle, Sh_cell{i,iM0}, iM0, NM, M0, dM);
        for k=1:NM
            Sh_cell{i,k}=Sh_tmp{k};
        end
    catch
    end
end


for i=fliplr(1:(iT0-1))
    fprintf('============================================================\n')
    fprintf('j = %d,  \t tau = %1.3f\n',i,TT(i)) %Print out progress
    try
    psimaxmin_in=[Sh_cell{i+1,iM0}.psimax,Sh_cell{i+1,iM0}.psimin];
    Sh_cell{i,iM0}=Sh_handle(Z,n,m, TT(i), M0, t, nu_star, psimaxmin_in, tol);
    Sh_tmp=Shock_pkg_new.Mach_sweep(Sh_handle, Sh_cell{i,iM0}, iM0, NM, M0, dM);
    for k=1:NM
        Sh_cell{i,k}=Sh_tmp{k};
    end
    catch
    end
end


end%end function

