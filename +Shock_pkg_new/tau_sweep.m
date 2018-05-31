function [Sh_cell] = tau_sweep(pre_calc_Shock, i0, N_steps, tau0, dtau)
% A function that takes a Shock function handle, and sweeps over tau
% to find the phimax and phimin as functions of Mach number.
% 
% The sweep is done by finding new shock solutions a step dtau away from a
% previously known shock solution. The function takes in a pre-calculated
% shock (pre_calc_Shock), the number of steps (N_steps), a starting tau
% (tau0), and a step size (dtau). It then runs up from i0+1 to N_steps
% followed by a run down from i0-1 to 1. 
% 
% If a problem is encountered on the way, dtau is halved until a solutions
% is found, when changing direction dM is resotred to its original value.


Sh_handle=@Shock_pkg_new.Shock_col;

%initializing the cell list to be returned
Sh_cell=cell(N_steps,1);
Sh_cell{i0}=pre_calc_Shock;

%Sweeping up
T_tmp=tau0;
dT_tmp=dtau;
for i=(i0+1):N_steps
    T_tmp=T_tmp+dT_tmp;%Stepping up the Mach #
    fprintf('============================================================\n')
    fprintf('i = %d,  \t tau = %1.3f\n',i,T_tmp) %Print out progress
    try
        [Sh_cell{i}, T_tmp, dT_tmp]=single_step(Sh_handle,Sh_cell{i-1}, T_tmp, dT_tmp, +1);
    catch
        warning('i = %d, failed. Now continuing with the down going.',i)
        break
    end
end

%Sweeping down
T_tmp=tau0;
dT_tmp=dtau;
for j=1:(i0-1)
    i=i0-j;
    T_tmp=T_tmp-dT_tmp;%Stepping down the Mach #
    fprintf('============================================================\n')
    fprintf('i = %d,  \t tau = %1.3f\n',i,T_tmp) %Print out progress
    Trying to find shock solution
    try
        [Sh_cell{i}, T_tmp, dT_tmp]=single_step(Sh_handle,Sh_cell{i+1}, T_tmp, dT_tmp, -1);
    catch
        warning('i = %d, failed. Returning.',i)
        return
    end
end 
fprintf('DONE!\n')%Extracting shock parameters
%tol=Sh_prev.tol;
%m=Sh_prev.m;
%n=Sh_prev.n;
%Z=Sh_prev.Z;
%tau=Sh_prev.tau;
%F_in=Sh_prev.F;


end%end function


function [Sh_tmp, T_tmp, dT_tmp] = single_step(Sh_handle,Sh_prev, T_tmp, dT_tmp, step_direction)
% A function for a single step.
% In this function the halving of dM, until successful, is implemented.

if isequal(Sh_handle,@Shock_tr)||isequal(Sh_handle,@Shock_pkg.Shock_tr)
    fprintf('Function not yet implemeted for trapping.\n\n')
    %args={Sh_prev.Z,Sh_prev.n,Sh_prev.tau,...
    %    Sh_prev.trapping_coef,T_tmp, [Sh_prev.psimax,Sh_prev.psimin], Sh_prev.tol};
    %i_T=4;
else
    args={Sh_prev.Z,Sh_prev.n,...
        T_tmp,Sh_prev.M, Sh_prev.t, Sh_prev.nu_star, [Sh_prev.psimax,Sh_prev.psimin], Sh_prev.tol};
    i_T=3;
end

%Trying to calculate a new shock
Sh_tmp=Sh_handle(args{:});


%If that fails try a tau closer to the previous (working) tau, by
%halving dtau
if (isnan(Sh_tmp.psimax) || isnan(Sh_tmp.psimin))&&refine_steps
    n_tries=0;
    while ( isnan(Sh_tmp.psimax) || isnan(Sh_tmp.psimin) )&& dT_tmp>=1e-7
        n_tries=n_tries+1;
        dT_tmp=dT_tmp/2; %decrease step size
        T_tmp=T_tmp-step_direction*dT_tmp; %Step back
        args{i_T}=T_tmp;
        fprintf('tau = %1.7f, dtau=%1.3e,   (n_tries=%d)\n',T_tmp,dT_tmp,n_tries)%Print out progress
        %Try to find shock solutions
        Sh_tmp=Sh_handle(args{:});

    end
    %decrease step size, to not re-try the previous, failed, Mach #
    dT_tmp=dT_tmp/2;
end %end if
end %end function