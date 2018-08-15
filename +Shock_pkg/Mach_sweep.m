function [Sh_cell] = Mach_sweep(Sh_handle, pre_calc_Shock, i0, N_steps, M0, dM, output_file)
% A function that takes a Shock function handle, and sweeps over Mach
% numbers to find the phimax and phimin as functions of Mach number.
% 
% The sweep is done by finding new shocksolutions a step dM away from a
% previously known shock solution. The function takes in a pre-calculated
% shock (pre_calc_Shock), the number of steps (N_steps), a sparting Mach
% number (M0), and a step size (dM). It then runs up from i0+1 to N_steps
% followed by a run down from i0-1 to 1.
% 
% If a problem is encountered on the way, dM is halved until a solutions is
% found, when changing direction dM is resotred to its original value.

if nargin<=6
    output_file='data/mach_sweep-tmp.tsv';
end
if exist(output_file, 'file')==2
    delete(output_file)
end


%initializing the cell list to be returned
Sh_cell=cell(N_steps,1);
Sh_cell{i0}=pre_calc_Shock;

data=[Sh_cell{i0}.Mach, Sh_cell{i0}.phimax, Sh_cell{i0}.phimin];
save(output_file, 'data', '-ascii','-append')

%Sweeping up
M_tmp=M0;
dM_tmp=dM;
for i=(i0+1):N_steps
    M_tmp=M_tmp+dM_tmp;%Stepping up the Mach #
    fprintf('============================================================\n')
    fprintf('i = %d,  \t Mach = %1.7f\n',i,M_tmp) %Print out progress
    try
        [Sh_cell{i}, M_tmp, dM_tmp]=single_step(Sh_handle,Sh_cell{i-1}, M_tmp, dM_tmp, +1);
        data=[Sh_cell{i}.Mach, Sh_cell{i}.phimax, Sh_cell{i}.phimin];
        save(output_file, 'data', '-ascii','-append')
    catch
        warning('i = %d, failed. Now continuing with the down going.',i)
        break
    end
end

%Sweeping down
M_tmp=M0;
dM_tmp=dM;
for j=1:(i0-1)
    i=i0-j;
    M_tmp=M_tmp-dM_tmp;%Stepping down the Mach #
    fprintf('============================================================\n')
    fprintf('i = %d,  \t Mach = %1.7f\n',i,M_tmp) %Print out progress
    %Trying to find shock solution
    try
        [Sh_cell{i}, M_tmp, dM_tmp]=single_step(Sh_handle,Sh_cell{i+1}, M_tmp, dM_tmp, -1);
        data=[Sh_cell{i}.Mach, Sh_cell{i}.phimax, Sh_cell{i}.phimin];
        save(output_file, 'data', '-ascii','-append')
    catch
        warning('i = %d, failed. Returning.',i)
        return
    end
end 
fprintf('DONE!\n')

end%end function


function [Sh_tmp, M_tmp, dM_tmp] = single_step(Sh_handle,Sh_prev, M_tmp, dM_tmp, step_direction)
% A function for a single step.
% In this function the halving of dM, until successful, is implemented.

%Extracting shock parameters
%tol=Sh_prev.tol;
%m=Sh_prev.m;
%n=Sh_prev.n;
%Z=Sh_prev.Z;
%tau=Sh_prev.tau;
%F_in=Sh_prev.F;


if isequal(Sh_handle,@Shock_tr)||isequal(Sh_handle,@Shock_pkg.Shock_tr)
    args={Sh_prev.m,Sh_prev.Z,Sh_prev.n,Sh_prev.tau,...
        Sh_prev.trapping_coef,'Mach',M_tmp, [Sh_prev.phimax,Sh_prev.phimin], Sh_prev.tol};
    i_M=7;
else
    args={Sh_prev.m,Sh_prev.Z,Sh_prev.n,Sh_prev.tau,...
        'Mach',M_tmp, [Sh_prev.phimax,Sh_prev.phimin], Sh_prev.tol};
    i_M=6;
end

%Trying to calculate a new shock
Sh_tmp=Sh_handle(args{:});
%Sh_tmp=Sh_handle(m,Z,n, tau,'Mach',M_tmp, F_in, tol);

%If that fails try a Mach # closer to the previous (working) Mach #, by
%halving dM
if repeat_cond(Sh_tmp) %isnan(Sh_tmp.phimax) || isnan(Sh_tmp.phimin) || Sh_tmp.Phi(-1,0)<0
    n_tries=0;
    if args{i_M+1}(2)<args{i_M+1}(1)*1e-2
        fprintf('Trying with a larger phi_maxmin_in.\nOld value = [%0.4f, %0.4f],\n',args{i_M+1}(1),args{i_M+1}(2))
        args{i_M+1}(2)=args{i_M+1}(1)*1e-2;
        fprintf('New value = [%0.4f, %0.4f].\n',args{i_M+1}(1),args{i_M+1}(2))
        Sh_tmp=Sh_handle(args{:});
        n_tries=1;
        if ~repeat_cond(Sh_tmp)
            return
        end
    end
    while repeat_cond(Sh_tmp) && dM_tmp>1e-5
        n_tries=n_tries+1;
        dM_tmp=dM_tmp/2; %decrease step size
        M_tmp=M_tmp-step_direction*dM_tmp; %Step back
        args{i_M}=M_tmp;
        fprintf('Mach = %1.7f, dM=%1.3e,   (n_tries=%d)\n',M_tmp,dM_tmp,n_tries)%Print out progress
        %Try to find shock solutions
        Sh_tmp=Sh_handle(args{:});
    end
    %decrease step size, to not re-try the previous, failed, Mach #
    dM_tmp=dM_tmp/2;
end %end if
end %end function

function R = repeat_cond(Sh)
R= isnan(Sh.phimax) || isnan(Sh.phimin) || Sh.Phi(-1,0)<0;
%The last condition (the one on Phi(-1,0)) is what causes the false error
%triggers for high tau sweeps.
end

