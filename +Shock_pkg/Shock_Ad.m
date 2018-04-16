classdef Shock_Ad < Shock_pkg.Shock_tr
% Shocks using the assumption that the trapped electrons are adabatic and
% therefore have flat distribution function in the trapped regions.
%
% 
%
% Author: Andréas Sundström
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There is still work that needs to be done on this code. For some reason
% we don't seem to be getting the right values of phimin all the time.

properties
        % This class has no own properties other than what is inhereted from
    % the Shock superclass.
end


methods
    function obj = Shock_Ad(m,Z,n, tau, speed_type,speed,phimaxmin_in, tol)
        %Constructor
        % Since there are no new properties, the constructor is mostly the same as
        % in the Shock super class
        import Shock_pkg.*
        if nargin==0
            args={};
        else
            args={m,Z,n, tau, [], speed_type,speed, phimaxmin_in,tol};
        end
        % Superclass construct
        obj=obj@Shock_pkg.Shock_tr(args{:});

        obj.trapping_coef=obj.phimin/obj.phimax;
    end %end Constructor
    
    function [n_el] = ne(obj, phi)
        % The total elctron density, due to ALL ion species.
        n_el=0; %init
        for i=1:obj.ion_species %summing over all species
        n_el=n_el+obj.ne_static(obj.m(i),obj.Z(i),obj.n(i),...
            phi, obj.phimax, obj.phimin, obj.V, obj.tau);
        end
    end

end %end methods

methods (Access=protected)
    function propgrp = getPropertyGroups(~)
        %Function defining how to display the properties
        proplist = {'tol','n','Z','m','ion_species','tau','Mach','V','cs',...
            'phimax','phimin', 'trapping_coef', 'F'};
        propgrp = matlab.mixin.util.PropertyGroup(proplist);
    end
    %{
    function [phimax, phimin]=find_phimaxmin(obj, F_in)
        % Finds the correct value of phimax, used in the constructor.
        
        % It is possible to specify starting points fr both phimax and
        % phimin, but if only one value is specified then 0.8 of the
        % initital guess for phimax is used for phimin. 
        if length(F_in)==2
            phiM_in=[F_in(1);F_in(2)]*(obj.Mach^2*obj.tau/2);
        else
            phiM_in=[1; 0.8]*F_in*(obj.Mach^2*obj.tau/2);
        end
        
        % The system of functions which we want to find the root of.
        %                   Phi_static(USDS,m,Z,n, phi, phimax, phimin, V, tau, tol)
        zerofun=@(phim)[obj.Phi_static(+1,obj.m,obj.Z,obj.n, phim(1), phim(1), phim(2), obj.V, obj.tau, obj.tol);
                        obj.Phi_static(-1,obj.m,obj.Z,obj.n, phim(2), phim(1), phim(2), obj.V, obj.tau, obj.tol)];
        % finding the root
        phiM=fsolve(zerofun, phiM_in, optimset('tolfun',obj.tol));
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
    %}
    function FVAL = zerofun_phiminmax(obj, phim)
        % The system of functions which we want to find the root of.
        % Phi_static(USDS,m,Z,n, phi, phimax, phimin, V, tau, tol)
        FVAL=[obj.Phi_static(+1,obj.m,obj.Z,obj.n, phim(1), phim(1), phim(2), obj.V, obj.tau, obj.tol);
              obj.Phi_static(-1,obj.m,obj.Z,obj.n, phim(2), phim(1), phim(2), obj.V, obj.tau, obj.tol)];
    end

end %end methods

%{
methods (Static=true, Access=protected)
    % Defining static functions which in turn are defined to be used in the
    % above functions, these should and can not be used out side this class
        
    function [ne_j] = ne_static(mj,Zj,nj, phi, phimax, phimin, V, tau)
        % The elctron density, due to ONE ion species. To get the total
        % value, sum this function over all ions. This is assuming that the
        % electron distribution function is flat in the trapped regions

        % The ion density (of this specific species) far upstream from the
        % shock, but still containing reflected ions. ("n" is without the
        % refelcted ions.) 
        njInfUS=0.5*nj*( 1 + 2*erf(V*sqrt(mj/2)) ...
            + erf((sqrt(2*Zj*phimax/mj)-V)*sqrt(mj/2)) );
        
        if phi<phimin
            ne_j=Zj*njInfUS*exp(phi/tau);
        else
            Y=real(sqrt((phi-phimin)/tau));
            ne_j=Zj*njInfUS*(2*Y*exp(phimin/tau)/sqrt(pi) + exp(phi/tau).*erfc(Y));
            % erfc(z)=1-erf(z)
        end
    end
    
    function [rho] = charge_dens_static(USDS, m,Z,n, phi, phimax, phimin, V, tau, tol)
        %  Note reqiures that m, Z, and n are all the same length.
        import Shock_pkg.Shock_Ad
        L=length(Z);
        rho=0;
        for j=1:L
        rho=rho+( -Shock_Ad.ne_static(m(j),Z(j),n(j), phi, phimax, phimin, V, tau)...
                  +Z(j)*Shock_Ad.nj_static(USDS,m(j),Z(j),n(j), phi, phimax, V, tol) );
        end
    end

    function [Phi_single] = Phi_static(USDS,m,Z,n, phi, phimax, phimin, V, tau, tol)
        %Must be static to be able to use this in find_phimax().
        import Shock_pkg.Shock_Ad
        % The integration limits are different for US and DS.
        if USDS==1
            Phi_single=integral(@(phiP) Shock_Ad.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, phimin, V, tau, tol), 0, phi, 'RelTol',tol);
        elseif USDS==-1
            Phi_single=integral(@(phiP) Shock_Ad.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, phimin, V, tau, tol), phimax, phi, 'RelTol',tol);
        else
            Phi_single=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
        end
    end
end %end methods
%}


end %end class