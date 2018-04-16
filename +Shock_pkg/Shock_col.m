classdef Shock_col < Shock_pkg.Shock_MB & Shock_pkg.Shock
% Shocks with a small collisionality. 
% 
% ONLY SINGLE ION-SPECIES PLASMAS!
%
% Author: Andréas Sundström
%

properties
    t
    nu_star
end

properties (Dependent)
    Upsilon
    phiA
end



methods
    function obj = Shock_col(m,Z,n, tau, speed_type,speed, t, nu_star, phimaxmin_in, tol)
        %Constructor
        % Since there are no new properties, the constructor is mostly the same as
        % in the Shock super class
        %import Shock_pkg.*
        if nargin==0
            args={};
        else
            args={m,Z,n, tau, speed_type,speed, tol};
        end
        % Superclass construct
        obj=obj@Shock_pkg.Shock(args{:});
        
        %{
        if length(phiminmax_prev)~=2
            warning('Bad phiminmax_prev.')
            return
        else
            obj.phimin_prev=min(phiminmax_prev);
            obj.phimax_prev=max(phiminmax_prev);
        end
        %}
        
        % If everything is ok, this finds and sets phimax
        if ( obj.ion_species>0 )&&( ~isempty(obj.V) )
            
            %if length(t)==1 || t(1)~=0
            %    t=[0,reshape(t,1,[])];
            %end
            obj.t=t; % At this stage length(t) must be 1
            obj.nu_star=nu_star;
            
            %Find phiminmax
            
            [obj.phimax,obj.phimin]=find_phimaxmin(obj, phimaxmin_in);
            
            % Checks that the solution is physical (actual maximum at
            % phi=phimax).
            if obj.charge_dens(+1, obj.phimax)<0 || obj.charge_dens(-1, obj.phimax)<0
                fprintf('ERROR: found phimax = %1.4f is NOT a true maximum.\n',obj.phimax)
                obj.phimax=NaN;
                return
            end
            
            %obj.phimin=find_phimin(obj, obj.phimin_prev);
            %if obj.phimin<0
                %obj.phimax=NaN;
            %    obj.phimin=NaN;
            %end
        end
    end %end Constructor
    
    function A=get.phiA(obj)
        % Get function for phiA (dependent variable)
        A=obj.phimax-obj.phimin;
    end
    function Ups=get.Upsilon(obj)
        % Get function for Upsilon (dependent variable)
        Ups=sqrt(2*obj.phiA/obj.nu_star);
        % Note, there is a calculation of Upsilon in zerofun_phiminmax(obj, phim)
    end
    
    function [ni] = nj(obj, USDS, phi, index)
        %Calculates the UpStream ion density of ion species number 'index'.
        %This is done using the static function nj_single.
        % The argument USDS must be either +1 (US) or -1 (DS).
        %fprintf('obj.nj\n')
        if abs(USDS)==1
            %obj.m(index),obj.Z(index),obj.n(index),
            %obj.phimax,obj.phimax_prev,obj.phimin_prev, obj.t,obj.Upsilon,
            %obj.V, USDS, obj.tol
                
            ni=arrayfun(@(x) obj.nj_single(obj.m(index),obj.Z(index),obj.n(index),...
                x, obj.phimax,obj.phimin, obj.t,obj.Upsilon,...
                obj.V, USDS, obj.tol), phi);            
        else
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            ni=NaN;
        end
    end
    %{
    function [ni_diff] = nj_diff(obj, USDS, phi, index)
        %Calculates the UpStream ion density of ion species number 'index'.
        %This is done using the static function nj_single.
        % The argument USDS must be either +1 (US) or -1 (DS).
        if abs(USDS)==1
            ni_diff=obj.nj_diff_static(USDS, obj.m(index),obj.Z(index),obj.n(index),...
                phi, obj.phimax,obj.phimax_prev,obj.phimin_prev, obj.t,obj.Upsilon, obj.V);
        else
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            ni_diff=NaN;
        end
    end
    %}
    function [n_el] = ne(obj, phi)
        % The total elctron density, due to ALL ion species.
        n_el=0; %init
        for i=1:obj.ion_species %summing over all species
        n_el=n_el+obj.ne_static(obj.m(i),obj.Z(i),obj.n(i),...
            phi, obj.phimax,obj.phimin, obj.t,obj.Upsilon, obj.V, obj.tau, obj.tol);
        end
    end

    
    function [X,phi,E,rho] = find_phi(obj,Xmin,Xmax)
        % This function calcualtes the potential phi as a function of x,
        % also the electric field E=dphi/dx comes "free of charge".
        
        % Functions specifying the derivatives of phi (G)
        odefunUS=@(x, G) [G(2);-obj.charge_dens(+1, G(1))];
        odefunDS=@(x, G) [G(2);-obj.charge_dens(-1, G(1))];
        %odefunUS=@(x, G) [G(2);-obj.charge_dens(+1, G(1))-obj.Z*obj.nj_diff(+1, G(1), 1)];
        %odefunDS=@(x, G) [G(2);-obj.charge_dens(-1, G(1))-obj.Z*obj.nj_diff(-1, G(1), 1)];
        G0=[obj.phimax;0]; %initial condition
        % If there is only one limit, then use a symmetrical inteval
        if nargin==2
            Xmax=-Xmin;
        end
        % Checks to see if the initial conditions are physical, i.e.
        % whether the electrostatic potential actually has a maximum at
        % phi=phimax.

        if ( sum(odefunUS(0,G0))<0 )&&( sum(odefunDS(0,G0))<0 )
            opt=odeset('reltol',obj.tol);
            [xUS, phiUS]=ode45(odefunUS, [0,Xmax],G0, opt);
            [xDS, phiDS]=ode45(odefunDS, [0,-abs(Xmin)],G0, opt);
            X  =[flip(xDS,1);xUS];
            phi=[flip(phiDS(:,1),1);phiUS(:,1)];
            E  =[flip(phiDS(:,2),1);phiUS(:,2)];
            fprintf('ODE solved!\n')
            
            %DEBUG
            rho=zeros(length(X),1);
            for i=1:length(rho)
                if X(i)>0
                    dG=odefunUS(X(i), [phi(i),E(i)]);
                else
                    dG=odefunDS(X(i), [phi(i),E(i)]);
                end
                rho(i)=-dG(2);
            end
        else
            X=[];phi=[];E=[];
            fprintf('ERROR: not a maximum at phi=phimax, charge density is negative.\n')
        end
    end
    
end %end methods

methods (Access=protected)
    function propgrp = getPropertyGroups(~)
        %Function defining how to display the properties
        proplist = {'tol','n','Z','m','ion_species','tau','Mach','V','cs',...
            'phimax','phimin','phiA','F', 'nu_star','t', 'Upsilon'};
        propgrp = matlab.mixin.util.PropertyGroup(proplist);
    end

    % {
    function [phimax, phimin]=find_phimaxmin(obj,phimaxmin_in)
        % Finds the correct value of phimax, used in the constructor.
        
        %{
        % It is possible to specify starting points for both phimax and
        % phimin, but if only one value is specified then 0.8 of the
        % initital guess for phimax is used for phimin. 
        if length(F_in)==2
            phiM_in=[F_in(1);F_in(2)]*(obj.Mach^2*obj.tau/2);
        else
            phiM_in=[1; 0.8]*F_in*(obj.Mach^2*obj.tau/2);
        end
        %}
        if length(phimaxmin_in)==2
            phiM_in=phimaxmin_in;
        else
            phiM_in=[1,0.7]*phimaxmin_in;
        end
        
        % finding the root
        phiM=fsolve(@(phim) obj.zerofun_phiminmax(phim), phiM_in, optimset('tolfun',obj.tol));
        % Error if there are some obviously wonky values
        if ( phiM(1)>phiM(2) )&& ( phiM(2)>0 )
            phimax=phiM(1);
            phimin=phiM(2);
        else
            phimax=NaN;
            phimin=NaN;
            fprintf('ERROR: something is wrong. phimax = %.02f, phimin = %0.2f.\n',...
                phiM(1),phiM(2))
        end
    end
    
    function FVAL = zerofun_phiminmax(obj, phim)
        % The system of functions which we want to find the root of.
        % Phi_static(USDS,m,Z,n, phi, phimax, phi_tr, V, tau, tol)
        Upsilon_guess=sqrt(2*(phim(1)-phim(2))/obj.nu_star);
        FVAL=[obj.Phi_static(+1,obj.m,obj.Z,obj.n, phim(1), phim(1), phim(2), obj.t,Upsilon_guess, obj.V, obj.tau, obj.tol);
              obj.Phi_static(-1,obj.m,obj.Z,obj.n, phim(2), phim(1), phim(2), obj.t,Upsilon_guess, obj.V, obj.tau, obj.tol)];
        % Note: using phim(1)*obj.trapping_coef, since at this stage,
        % phimax has not yet been set.
    end
    % }
    %{
    function [phimax]=find_phimax(obj)
        % Finds the correct value of phimax, used in the constructor.
        phimaxin=obj.phimax_prev;
        %The Sagdeev potential has to be 0 at phi=phimax
        
        phimax=fzero(@(phim) obj.Phi_static(+1,obj.m,obj.Z,obj.n,...
                        phim, phim, obj.phimax_prev,obj.phimin_prev, ...
                        obj.t,obj.Upsilon, obj.V, obj.tau, obj.tol), phimaxin,...
                        optimset('tolfun',obj.tol));
    end
    
    function [phimin]=find_phimin(obj, phimin_in)
        % Finds phimin is obj is a valid shock.
        %Checks to see if all is well
        
        %fprintf('find_phimin\n') %DEBUG
        if isnan(obj.phimax)
            phimin=NaN;
            return
        end
        
        %The DS Sagdeev potential has to be 0 at phi=phimmin
        phimin=fzero(@(phim) obj.Phi_static(-1,obj.m,obj.Z,obj.n,...
                    phim, obj.phimax, obj.phimax_prev,obj.phimin_prev, ...
                    obj.t,obj.Upsilon, obj.V, obj.tau, obj.tol), phimin_in,...
                    optimset('tolfun',obj.tol));
        if phimin>=obj.phimax
            fprintf('The specified phimin_in found the wrong root of Phi.\nNow trying in the interval [0, 1-tol]*phimax.\n\n')
            %Now we instead define an interval in which to look
            phimin_in=[0,(1-obj.tol)]*obj.phimax;
            if obj.Phi(-1, phimin_in )*obj.Phi(-1, 0)<0
                phimin=find_phimin(obj, phimin_in);
            else
                phimin=NaN;
                fprintf('ERROR: No phimin in the interval [0, 1-tol]*phimax.\n')
            end
        end
    end
    %}
    
    function [PHI] = Phi_single(obj, USDS, phi)
        %Functino that takes a single phi value, used in Phi().
        if USDS==1
            PHI=integral(@(phiP) obj.charge_dens(USDS, phiP), 0, phi, 'RelTol',obj.tol);
        elseif USDS==-1
            PHI=integral(@(phiP) obj.charge_dens(USDS, phiP), obj.phimax, phi, 'RelTol',obj.tol);
            %PHI=PHI-obj.Phi_diff_static(USDS, obj.m,obj.Z,obj.n, obj.phimax, obj.phimax,...
            %obj.phimax_prev,obj.phimin_prev, obj.t,obj.Upsilon, obj.V);
        else
            PHI=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            return
        end
        PHI=real(PHI);
        %+obj.Phi_diff_static(USDS, obj.m,obj.Z,obj.n, phi, obj.phimax,...
        %    obj.phimax_prev,obj.phimin_prev, obj.t,obj.Upsilon, obj.V);
    end
    
end %end methods (proteccted)

methods (Static=true, Access=protected)
    % Defining static functions which in turn are defined to be used in the
    % above functions, these should and can not be used out side this class
    
    
    function [Phi_single] = Phi_static(USDS,m,Z,n, phi, phimax, phimin, t,Upsilon, V, tau, tol)
        %Must be static to be able to use this in find_phimax().
        import Shock_pkg.Shock_col
        if USDS==1
            Phi_single=integral(@(phiP) Shock_col.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, phimin,...
                t,Upsilon, V, tau, tol), 0, phi, 'RelTol',tol);
        elseif USDS==-1
            Phi_single=integral(@(phiP) Shock_col.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, phimin,...
                t,Upsilon, V, tau, tol), phimax, phi, 'RelTol',tol);
        else
            Phi_single=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            return
        end
        
        Phi_single=real(Phi_single);
    end

    
    function [rho] = charge_dens_static(USDS, m,Z,n, phi, phimax,phimin, t,Upsilon, V, tau, tol)
        %  Note reqiures that m, Z, and n are all the same length.
        import Shock_pkg.Shock_col
        L=length(Z);
        rho=0;
        %fprintf('rho\n') %DEBUG
        for j=1:L
        rho=rho+( +Z(j)*Shock_col.nj_static(USDS,m(j),Z(j),n(j), phi, phimax,phimin, t,Upsilon, V, tol)...
                            -Shock_col.ne_static(m(j),Z(j),n(j), phi, phimax,phimin, t,Upsilon, V, tau, tol));
        end
    end
    
    %%{
    function [n] = nj_static(USDS,mj,Zj,nj, phi, phimax,phimin, t,Upsilon, V, tol)
        %Calclulats the upstream or downstream ion density for ONE ion
        %species, due to the value of the potential and the maximum value
        %of the potential. 
        import Shock_pkg.Shock_col
        n=arrayfun(@(x) Shock_col.nj_single(mj,Zj,nj, x, phimax,phimin, t,Upsilon, V, USDS, tol), phi);
    end
    
    function [n] = nj_single(mj,Zj,nj, phi, phimax, phimin, t,Upsilon, V, USDS, tol)
        %Helping function to nj_static
        %This function is needed becase integral() passes a phi as a vector
        %argument, while in turn integral() itselv can't handel vectors as
        %endpoints.
        % the varible lim is either +1 or -1.
        import Shock_pkg.Shock_col
        v0=real(sqrt(2*Zj*(phimax-phi)/mj));
        n=integral(@(v) Shock_col.fj_static(mj,Zj,nj, phi, V, v), -Inf,USDS*v0, 'RelTol', tol);
        
        
        %Here is where the collision come in
        phiA=abs(phimax-phimin);
        %n_reg=2-USDS;
        lim=10*sqrt(t)/Upsilon;
        int_diff=integral(@(k) Shock_col.fjIII_reduced(k, phi, phimax,phimin, t,Upsilon),...
             0,sqrt(2*lim+lim^2), 'RelTol', tol);
             %1,1+lim, 'RelTol', tol);
        
        if USDS == -1
            Qphi=sqrt((phimax-phi)/phiA);
            int_diff=int_diff+ 2*integral(@(k) Shock_col.fjII_reduced(k, phi, phimax,phimin, t,Upsilon),...
             0,Qphi, 'RelTol', tol);
            
        end
        n=n+ real(int_diff).*Shock_col.fj_static(mj,Zj,nj, phi, V, -v0)*sqrt(2*phiA);
    end
    
    function f3 = fjIII_reduced(k, phi, phimax,phimin, t,Upsilon)
        srt2=.5/sqrt(t);
        %{
        phiA=abs(phimax-phimin_prev);
        phitilde=phi-phimin_prev;

        f3=erfc(-Shock_pkg.w_approx(1-k,Upsilon)*srt2).*real(k./sqrt(k.^2-phitilde/phiA));
        %}
        % {
        phiA=abs(phimax-phimin);
        Qphisq=(phimax-phi)/phiA;
        if Qphisq~=0
            f3=erfc(-Shock_pkg.w_approx(1-real(sqrt(k.^2+1)),Upsilon)*srt2).*real(k./sqrt(k.^2+Qphisq));
        else
            f3=erfc(-Shock_pkg.w_approx(1-real(sqrt(k.^2+1)),Upsilon)*srt2);
        end
        % }
    end
    function f2 = fjII_reduced(k, phi, phimax,phimin, t,Upsilon)
        srt2=.5/sqrt(t);
        phiA=abs(phimax-phimin);
        Qphisq=(phimax-phi)/phiA;
        if Qphisq>0
            f2=erfc(Shock_pkg.w_approx(1-real(sqrt(1-k.^2)),Upsilon)*srt2).*real(k./sqrt(Qphisq-k.^2));
        else
            f2=zeros(size(k));
        end
        % }
    end
    %%}
    
    function [ne_j] = ne_static(mj,Zj,nj, phi, phimax, phimin, t,Upsilon, V, tau, tol)
        % The elctron density, due to ONE ion species. To get the total
        % value, sum this function over all ions. This is assuming that the
        % electron distribution function is walways a pure
        % Maxwell-Boltzmann. 
        % The ion density (of this specific species) far upstream from the
        % shock, but still containing reflected ions. ("n" is without the
        % refelcted ions.) 
        %njInfUS=0.5*nj*( 1 + 2*erf(V*sqrt(mj/2)) +...
        %    erf(real(sqrt(2*Zj*phimax/mj)-V)*sqrt(mj/2)) );
        
        import Shock_pkg.Shock_col

        njInfUS=Shock_col.nj_static(+1,mj,Zj,nj,...
            0,phimax,phimin, t,Upsilon, V, tol);
        %njInfUS=njInfUS + Shock_col.nj_diff_static(+1,mj,Zj,nj,...
        %    0,phimax,phimax_prev,phimin_prev, t,Upsilon, V);
        
        % Given a pure M-B, the distribution integrates up to this density:
        ne_j=Zj*njInfUS*exp(phi/tau);
    end

end %end methods



end %end class