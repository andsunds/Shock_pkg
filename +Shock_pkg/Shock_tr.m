classdef Shock_tr < Shock_pkg.Shock
% Shocks using the assumption that the electrons get trapped at a level
% trapping_coef*phimax and therefore have flat distribution function in the
% trapped regions. 
%
% 
%
% Author: Andréas Sundström (c)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There is still work that needs to be done on this code. For some reason
% we don't seem to be getting the right values of phimin all the time.

properties (SetAccess = protected)
    trapping_coef
end
properties (Dependent)
    phi_tr
end

methods
    function obj = Shock_tr(m,Z,n,tau, C_tr, speed_type,speed,phimaxmin_in, tol)
        %Constructor
        import Shock_pkg.*
        if nargin==0
            args={};
        else
            args={m,Z,n, tau, speed_type,speed, tol};
        end

        % Superclass construct
        obj=obj@Shock_pkg.Shock(args{:});
        if nargin==9
            obj.trapping_coef=C_tr;
            % If everything is ok, this finds and sets phimax
            if ( obj.ion_species>0 )&&( ~isempty(obj.V) )
                [obj.phimax, obj.phimin]=find_phimaxmin(obj, phimaxmin_in);
            end
        else
            return
        end
        if obj.charge_dens(+1, obj.phimax)<0 || obj.charge_dens(-1, obj.phimax)<0
            fprintf('ERROR: found phimax = %1.4f is NOT a true maximum.\n',obj.phimax)
            obj.phimax=NaN;
            obj.phimin=NaN;
            return
        end
    end %end Constructor
    
    function phitr=get.phi_tr(obj)
        % Get function for the speed of sound (dependent variable)
        phitr=obj.trapping_coef*obj.phimax;
    end
    
    function [n_el] = ne(obj, phi)
        % The total elctron density, due to ALL ion species.
        n_el=0; %init
        for i=1:obj.ion_species %summing over all species
        n_el=n_el+obj.ne_static(obj.m(i),obj.Z(i),obj.n(i),...
            phi, obj.phimax, obj.phi_tr, obj.V, obj.tau);
        end
    end

end %end methods

methods (Access=protected)
    function propgrp = getPropertyGroups(~)
        %Function defining how to display the properties
        proplist = {'tol','n','Z','m','tau','Mach','V','cs',...
            'phimax','phimin','phi_tr', 'trapping_coef'};
        propgrp = matlab.mixin.util.PropertyGroup(proplist);
    end
   
    function [phimax, phimin]=find_phimaxmin(obj, phimaxmin_in)
        % Finds the correct value of phimax, used in the constructor.
        
        % It is possible to specify starting points for both phimax and
        % phimin, but if only one value is specified then 0.8 of the
        % initital guess for phimax is used for phimin. 
        if length(phimaxmin_in)==2
            phiM_in=[phimaxmin_in(1);phimaxmin_in(2)];%*(obj.Mach^2*obj.tau/2);
        else
            phiM_in=[1; 0.8]*phimaxmin_in;%*(obj.Mach^2*obj.tau/2);
        end
        
        % finding the root
        OPT=optimset('tolfun',obj.tol,'Display','off');
        [phiM,~,EF]=fsolve(@(phim) obj.zerofun_phiminmax(phim), phiM_in, OPT);
        fprintf('EF = %d.\n',EF)
        % Error if there are some obviously wonky values
        if  phiM(1)>phiM(2) && phiM(2)>0 
            phimax=phiM(1);
            phimin=phiM(2);
        else
            phimax=NaN;
            phimin=NaN;
            fprintf('ERROR: something is wrong. phimax = %.02f, phimin = %0.2f.\n',...
                phiM(1),phiM(2))
        end
        if EF<=0
            fprintf('No solution found, ExitFlag %d.\n',EF)
            phimax=NaN;
            phimin=NaN;
        end
    end
    
    function FVAL = zerofun_phiminmax(obj, phim)
        % The system of functions which we want to find the root of.
        % Phi_static(USDS,m,Z,n, phi, phimax, phi_tr, V, tau, tol)
        FVAL=[obj.Phi_static(+1,obj.m,obj.Z,obj.n, phim(1), phim(1), phim(1)*obj.trapping_coef, obj.V, obj.tau, obj.tol);
              obj.Phi_static(-1,obj.m,obj.Z,obj.n, phim(2), phim(1), phim(1)*obj.trapping_coef, obj.V, obj.tau, obj.tol)];
        % Note: using phim(1)*obj.trapping_coef, since at this stage,
        % phimax has not yet been set.
    end
    
end %end methods


methods (Static=true, Access=protected)
    % Defining static functions which in turn are defined to be used in the
    % above functions, these should and can not be used out side this class
        
    function [ne_j] = ne_static(mj,Zj,nj, phi, phimax, phi_tr, V, tau)
        % The elctron density, due to ONE ion species. To get the total
        % value, sum this function over all ions. This is assuming that the
        % electron distribution function is flat in the trapped regions

        % The ion density (of this specific species) far upstream from the
        % shock, but still containing reflected ions. ("n" is without the
        % refelcted ions.) 
        njInfUS=0.5*nj*( 1 + 2*erf(V*sqrt(mj/2)) ...
            + erf(real(sqrt(2*Zj*phimax/mj)-V)*sqrt(mj/2)) );
        
        if phi<phi_tr
            ne_j=Zj*njInfUS*exp(phi/tau);
        else
            Y=real(sqrt((phi-phi_tr)/tau));
            ne_j=Zj*njInfUS*(2*Y*exp(phi_tr/tau)/sqrt(pi) + exp(phi/tau).*erfc(Y));
            % erfc(z)=1-erf(z)
        end
    end
    
    function [rho] = charge_dens_static(USDS, m,Z,n, phi, phimax, phi_tr, V, tau, tol)
        %  Note reqiures that m, Z, and n are all the same length.
        import Shock_pkg.Shock_tr
        L=length(Z);
        rho=0;
        for j=1:L
        rho=rho+( -Shock_tr.ne_static(m(j),Z(j),n(j), phi, phimax, phi_tr, V, tau)...
                  +Z(j)*Shock_tr.nj_static(USDS,m(j),Z(j),n(j), phi, phimax, V, tol) );
        end
    end

    function [Phi_single] = Phi_static(USDS,m,Z,n, phi, phimax, phi_tr, V, tau, tol)
        %Must be static to be able to use this in find_phimax().
        import Shock_pkg.Shock_tr
        % The integration limits are different for US and DS.
        if USDS==1
            Phi_single=integral(@(phiP) Shock_tr.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, phi_tr, V, tau, tol), 0, phi, 'RelTol',tol);
        elseif USDS==-1
            Phi_single=integral(@(phiP) Shock_tr.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, phi_tr, V, tau, tol), phimax, phi, 'RelTol',tol);
        else
            Phi_single=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
        end
    end
end %end methods



end %end class