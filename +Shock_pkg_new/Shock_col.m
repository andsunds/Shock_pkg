classdef Shock_col < Shock_pkg_new.Shock
% Shocks with a small collisionality. 
% 
% ONLY SINGLE ION-SPECIES PLASMAS!
%
% Author: Andréas Sundström (c)
%

properties
    t
    nu_star
end

properties (Dependent)
    Upsilon
    psiA
end



methods
    function obj = Shock_col(Z,n, taui, Mach, t, nu_star, psimaxmin_in, tol)
        %Constructor
        % Since there are no new properties, the constructor is mostly the same as
        % in the Shock super class
        %import Shock_pkg.*
        if nargin==0
            args={};
        else
            args={Z,n, taui, Mach, tol};
        end
        
        % Superclass construct
        obj=obj@Shock_pkg_new.Shock(args{:});
        
        
        % If everything is ok, this finds and sets phimax
        if ( ~isempty(obj.M) )
            obj.t=t; 
            obj.nu_star=nu_star;
            
            %Find phiminmax
            [obj.psimax,obj.psimin]=find_psimaxmin(obj, psimaxmin_in);
            
            % Checks that the solution is physical (actual maximum at
            % phi=phimax).
            if obj.charge_dens(+1, obj.psimax)<0 || obj.charge_dens(-1, obj.psimax)<0
                fprintf('ERROR: found psimax = %1.4f is NOT a true maximum.\n',obj.psimax)
                obj.psimax=NaN;
                obj.psimin=NaN;
                return
            end
        end
    end %end Constructor
    
    function A=get.psiA(obj)
        % Get function for phiA (dependent variable)
        A=obj.psimax-obj.psimin;
    end
    function Ups=get.Upsilon(obj)
        % Get function for Upsilon (dependent variable)
        Ups=sqrt(2*obj.psiA*obj.taui/obj.nu_star);
        % Note, there is a calculation of Upsilon in zerofun_phiminmax(obj, phim)
    end
    

    function [ni] = nj(obj, USDS, psi, index)
        %Calculates the UpStream ion density of ion species number 'index'.
        %This is done using the static function nj_single.
        % The argument USDS must be either +1 (US) or -1 (DS).
        if nargin <4
            index=1;
        end
        if abs(USDS)==1 
            ni=arrayfun(@(x) obj.nj_single(obj.taui(index),obj.n(index),...
                x, obj.psimax,obj.psimin, obj.t,obj.Upsilon,...
                obj.M, USDS, obj.tol), psi);            
        else
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            ni=NaN;
        end
    end
    
    function [n_el] = ne(obj, psi)
        % The total elctron density, due to ALL ion species.
        n_el=0; %init
        for i=1:length(obj.n) %summing over all species
            n_el=n_el+obj.ne_static(obj.taui,obj.n(i),...
                psi, obj.psimax,obj.psimin, obj.t,obj.Upsilon, obj.M, obj.tol);
        end
    end

 
end %end methods

methods (Access=protected)
    function propgrp = getPropertyGroups(~)
        %Function defining how to display the properties
        proplist = {'tol','n','Z','m','ion_species','taui','M',...
            'psimax','psimin','psiA', 'nu_star','t', 'Upsilon'};
        propgrp = matlab.mixin.util.PropertyGroup(proplist);
    end

    % {
    function [psimax, psimin]=find_psimaxmin(obj,psimaxmin_in)
        % Finds the correct value of phimax, used in the constructor.
        
        %Depending on if the initial guess is for both phimax and phimin or
        %just phimax
        if length(psimaxmin_in)==2
            psiM_in=psimaxmin_in;
        else
            psiM_in=[1,0.7]*psimaxmin_in;
        end
        
        % finding the root
        psiM=fsolve(@(psim) obj.zerofun_psiminmax(psim), psiM_in, optimset('tolfun',obj.tol));
        % Error if there are some obviously wonky values
        if ( psiM(1)>psiM(2) )&& ( psiM(2)>0 )
            psimax=psiM(1);
            psimin=psiM(2);
        else
            psimax=NaN;
            psimin=NaN;
            fprintf('ERROR: something is wrong. phimax = %.02f, phimin = %0.2f.\n',...
                psiM(1),psiM(2))
        end
    end
    
    function FVAL = zerofun_psiminmax(obj, psim)
        % The system of functions which we want to find the root of.
        % Phi_static(USDS,m,Z,n, phi, phimax, phi_tr, V, tau, tol)
        Upsilon_guess=sqrt(2*(psim(1)-psim(2))*obj.taui/obj.nu_star);
        FVAL=[obj.Psi_static(+1,obj.Z,obj.n, psim(1), psim(1), psim(2), obj.t,Upsilon_guess, obj.M, obj.taui, obj.tol);
              obj.Psi_static(-1,obj.Z,obj.n, psim(2), psim(1), psim(2), obj.t,Upsilon_guess, obj.M, obj.taui, obj.tol)];
        % Note: using phim(1)*obj.trapping_coef, since at this stage,
        % phimax has not yet been set.
    end

    
    function [PSI] = Psi_single(obj, USDS, psi)
        %Functino that takes a single phi value, used in Phi().
        if USDS==1 %Upstream
            PSI=integral(@(psiP) obj.charge_dens(USDS, psiP), 0, psi, 'RelTol',obj.tol);
        elseif USDS==-1 %Downstream
            PSI=integral(@(psiP) obj.charge_dens(USDS, psiP), obj.psimax, psi, 'RelTol',obj.tol);
        else
            PSI=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            return
        end
        PSI=real(PSI);
    end
end %end methods (proteccted)

methods (Static=true, Access=protected)
    % Defining static functions which in turn are defined to be used in the
    % above functions, these should and can not be used out side this class
    
    
    function [Psi_single] = Psi_static(USDS,Z,n, psi, psimax, psimin, t,Upsilon, M, taui, tol)
        %Must be static to be able to use this in find_phimax().
        import Shock_pkg_new.Shock_col
        if USDS==1 %Upstream
            Psi_single=integral(@(psiP) Shock_col.charge_dens_static(...
                USDS,Z,n, psiP, psimax, psimin,...
                t,Upsilon, M, taui, tol), 0, psi, 'RelTol',tol);
        elseif USDS==-1 %Downstream
            Psi_single=integral(@(psiP) Shock_col.charge_dens_static(...
                USDS,Z,n, psiP, psimax, psimin,...
                t,Upsilon, M, taui, tol), psimax, psi, 'RelTol',tol);
        else
            Psi_single=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            return
        end
        % Removes unphyical imaginary part of the numrical calculation
        Psi_single=real(Psi_single);
    end

    
    function [rho] = charge_dens_static(USDS, Z,n, psi, psimax,psimin, t,Upsilon, M, taui, tol)
        %  Note reqiures that m, Z, and n are all the same length.
        import Shock_pkg_new.Shock_col
        L=length(Z);
        rho=0;
        for j=1:L
        rho=rho+( Z(j)*Shock_col.nj_static(USDS,taui,n(j), psi, psimax,psimin, t,Upsilon, M, tol) ...
                 -Shock_col.ne_static(taui,n(j), psi, psimax,psimin, t,Upsilon, M, tol));
        end
        %rho=rho+( +Z(j)*Shock_col.nj_static(USDS,m(j),Z(j),n(j), phi, phimax,phimin, t,Upsilon, V, tol)...
        %                    -Shock_col.ne_static(m(j),Z(j),n(j), phi, phimax,phimin, t,Upsilon, V, tau, tol));
    end
    
    %%{
    function [n] = nj_static(USDS,taui,nj, psi, psimax,psimin, t,Upsilon, M, tol)
        %Calclulats the upstream or downstream ion density for ONE ion
        %species, due to the value of the potential and the maximum value
        %of the potential. 
        import Shock_pkg_new.Shock_col
        n=arrayfun(@(x) Shock_col.nj_single(taui,nj, x, psimax,psimin, t,Upsilon, M, USDS, tol),psi);
        %(mj,Zj,nj, x, phimax,phimin, t,Upsilon, V, USDS, tol), phi);
    end

    function [n] = nj_single(taui,nj, psi, psimax,psimin, t,Upsilon, M, USDS, tol)
        %Helping function to nj_static
        %This function is needed becase integral() passes a phi as a vector
        %argument, while in turn integral() itselv can't handel vectors as
        %endpoints.
        % the varible lim is either +1 or -1.
        import Shock_pkg_new.Shock_col
        v0=real(sqrt(2*(psimax-psi)));
        n=integral(@(v) Shock_col.fj_static(taui,nj, psi, M, v), -Inf,USDS*v0, 'RelTol', tol);
        
        if t>0
            %Here is where the contribution from the collision are added
            psiA=abs(psimax-psimin);%Amplitude of the DS oscillation
            %lim1=10*sqrt(t)/Upsilon;
            %lim=sqrt(2*lim1+lim1^2); %integration limit to help the num integration
            
            int_diff=integral(@(k) Shock_col.fjIII_reduced(k, psi, psimax,psimin, t,Upsilon),...
                 0,inf, 'RelTol', tol);
            
            if USDS == -1
                %In the downstream, we also have to include the trapped ions.
                Qphi=sqrt((psimax-psi)/psiA);
                int_diff=int_diff+ 2*integral(@(k) Shock_col.fjII_reduced(k, psi, psimax,psimin, t,Upsilon),...
                                                0,Qphi, 'RelTol', tol);
            end
            n=n+ real(int_diff)*sqrt(2*psiA).*Shock_col.fj_static(taui,nj, psi,M,-v0);%(mj,Zj,nj, phi, V, -v0)
        end
    end

    function f3 = fjIII_reduced(k, psi, psimax,psimin, t,Upsilon)
        %The ion distribution function in region III, co-passing region
        %(overtaking the shock). This distribution is present in both the
        %up- and the downstream.
        srt2=.5/sqrt(t);
        psiA=abs(psimax-psimin);
        Qphisq=(psimax-psi)/psiA;
        if Qphisq~=0
            %Shock_pkg.w_approx(1-real(sqrt(k.^2+1)),Upsilon)
            f3=erfc(-Shock_pkg.w_approx(1-real(sqrt(k.^2+1)),Upsilon)*srt2).*real(k./sqrt(k.^2+Qphisq));
        else
            f3=erfc(-Shock_pkg.w_approx(1-real(sqrt(k.^2+1)),Upsilon)*srt2);
        end
    end
    function f2 = fjII_reduced(k, psi, psimax,psimin, t,Upsilon)
        %The ion distribution function in region II, trapped region. This
        %population is only present in the downstream.
        srt2=.5/sqrt(t);
        psiA=abs(psimax-psimin);
        Qphisq=(psimax-psi)/psiA;
        if Qphisq>0
            f2=erfc(Shock_pkg.w_approx(1-real(sqrt(1-k.^2)),Upsilon)*srt2).*real(k./sqrt(Qphisq-k.^2));
        else
            f2=zeros(size(k));
        end
    end
    


    function [ne_j] = ne_static(taui,nj, psi, psimax, psimin, t,Upsilon, M, tol)
        % The elctron density, due to ONE ion species. To get the total
        % value, sum this function over all ions. This is assuming that the
        % electron distribution function is walways a pure
        % Maxwell-Boltzmann.
        import Shock_pkg_new.Shock_col

        %The ion density in the far upstream
        njInfUS=Shock_col.nj_static(+1,taui,nj, 0, psimax,psimin, t,Upsilon, M, tol);
        %nj_static(+1,mj,Zj,nj, 0,phimax,phimin, t,Upsilon, V, tol);
        
        % Given a pure M-B, the distribution integrates up to this density:
        ne_j=njInfUS*exp(psi);
    end

end %end methods



end %end class