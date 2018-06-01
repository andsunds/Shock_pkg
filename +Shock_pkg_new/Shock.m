classdef (Abstract=true) Shock < matlab.mixin.CustomDisplay
% Base class for Shocks. This is an abstrac class.
% This class contains some of the basic properties and methods needed to
% represent the 1D electrostatic shocks in plasmas studied in this project.
%
% The inheritance from matlab.mixin.CustomDisplay, is just so that the
% property display can be changed.
%
% Author: Andréas Sundström
%
% Superclass for: Shock_MB
    
properties (SetAccess = protected)
    tol  %The tolerance to which the integrations should be carried
    
    %These properties can be vectors, but they have to be the same length.
    Z
    n
    m
    taui         % electron-ion temp ration
    M           % flow speed
    psimax      % Max value of the electrostatic potential
    psimin % Min value of the electrostatic potential in downstream region
    
end
properties (Dependent)
    %lambda %Wavelength of DS oscillations
    N_ion % the number of ion species
    zeta
end
    
methods
    function obj = Shock(Z,n,m, taui, Mach, tol)
        %Constructor
        if nargin==0 %this is due to the way you call a superclass constructor
            %Do nothing
        else
            %Check that n and Z have all ion-species 
            if length(Z)==length(n) && length(Z)==length(m) && length(Z)==length(taui)
                %Set values upon construction.
                obj.tol=tol;
                obj.taui=taui;  %electron-ion temp ration
                obj.n=n/sum(n.*Z);  %ion densities (normalized so that main ion density is 1)
                obj.Z=Z;       %charges
                obj.m=m;
                obj.M=Mach;
            else
                warning('Z, n, m, and taui not the same length!')
            end
        end
    end %end consructor
    
    function Nion=get.N_ion(obj)
        Nion=length(obj.Z);
    end
    
    function z=get.zeta(obj)
        z=obj.Z*obj.m(1)./(obj.m*obj.Z(1));
    end
    
    function L=lambda(obj)
        %Get function for lambda
        L = sqrt(2)*integral(@(p)1./sqrt(-obj.Psi(-1,p)),...
            obj.psimin, obj.psimax, 'RelTol', max([obj.tol,1e-8]));
    end
       

    function [X,phi,E,rho] = find_psi(obj,Xmin,Xmax)
        % This function calcualtes the potential phi as a function of x,
        % also the electric field E=dphi/dx comes "free of charge".
        
        % Functions specifying the derivatives of phi (G)
        odefunUS=@(x, G) [G(2);-obj.charge_dens(+1, G(1))];
        odefunDS=@(x, G) [G(2);-obj.charge_dens(-1, G(1))];
        G0=[obj.psimax;0]; %initial condition
        % If there is only one limit, then use a symmetrical inteval
        if nargin==2
            Xmax=-Xmin;
        end
        % Checks to see if the initial conditions are physical, i.e.
        % whether the electrostatic potential actually has a maximum at
        % phi=phimax.
        if ( sum(odefunUS(0,G0))<0 )&&( sum(odefunDS(0,G0))<0 )
            opt=odeset('reltol',obj.tol);
            [xUS, psiUS]=ode45(odefunUS, [0,Xmax],G0, opt);
            [xDS, psiDS]=ode45(odefunDS, [0,-abs(Xmin)],G0, opt);
            X  =[flip(xDS,1);xUS];
            phi=[flip(psiDS(:,1),1);psiUS(:,1)];
            E  =[flip(psiDS(:,2),1);psiUS(:,2)];
            
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
    

    function [ni] = nj(obj, USDS, psi, index)
        %Calculates the UpStream ion density of ion species number 'index'.
        %This is done using the static function nj_single.
        % The argument USDS must be either +1 (US) or -1 (DS).
        if nargin <4
            index=1;
        end
        if abs(USDS)==1
            ni=arrayfun(@(x) obj.nj_single(obj.taui(index),obj.n(index),obj.zeta(index),...
                x, obj.psimax, obj.M, USDS, obj.tol), psi);
        else
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            ni=NaN;
        end
    end
    

    function [rho] = charge_dens(obj, USDS, psi)
        % Retruns the charge density either US or DS at potential phi .
        % The argument USDS must be either +1 (US) or -1 (DS).
        if abs(USDS)==1
            rho=obj.charge_dens_static(USDS, obj.Z,obj.n,obj.zeta, psi, obj.psimax, obj.M, obj.taui, obj.tol);
        else 
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            rho=NaN;
        end
    end

    function [PSI] = Psi(obj, USDS, psi)
        % Calculates the Sagdeev potential.
        % As per the theory, the Sagdeev potential is just the chagre
        % density integrated over the electrostatic potential. This is used
        % to calculate the maximum value of phi. 
        %    The argument USDS must be either +1 (US) or -1 (DS).
        if abs(USDS)==1 
            PSI=arrayfun(@(x) obj.Psi_single(USDS,x),psi);
        else
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            PSI=NaN;
        end
    end
end  %end methods

methods (Access=protected)
    
    function [PSI] = Psi_single(obj, USDS, psi)
        %Functino that takes a single phi value, used in Phi().
        if USDS==1
            PSI=integral(@(phiP) obj.charge_dens(USDS, phiP), 0, psi, 'RelTol',obj.tol);
        elseif USDS==-1
            PSI=integral(@(phiP) obj.charge_dens(USDS, phiP), obj.psimax, psi, 'RelTol',obj.tol);
        else
            PSI=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
        end
    end
end



methods (Abstract=true, Static=true, Access=protected)
    [ne_j] = ne_static(taui, Zj,nj,zetaj, psi, psimax, M);
    [rho] = charge_dens_static(USDS, Z,n,zeta, psi, psimax, M, taui, tol);
    [Psi_single] = Psi_static(USDS,Z,n,zeta, psi, psimax, M, taui, tol);
end



methods (Static=true, Access=protected)
    % Defining static functions which in turn are used in the definitions
    % of some of the above functions, these should and can not be used
    % outside of this class or its subclasses. 
    
    function [dist] = fj_static(taui,nj,zetaj, psi, M, v)
        % The ion distributin function of a specific ion species
        dist = nj*sqrt(taui/(2*pi)) * exp(-(taui/2)*(sqrt(v.^2+2*zetaj*psi)-M).^2);
    end
    function [n] = nj_static(USDS,taui,nj,zetaj, psi, psimax, M, tol)
        %Calclulats the upstream or downstream ion density for ONE ion
        %species, due to the value of the potential and the maximum value
        %of the potential. 
        import Shock_pkg_new.Shock
        %taui,nj, psi, psimax, M, USDS, tol
        n=arrayfun(@(x) Shock.nj_single(taui,nj,zetaj, x, psimax, M, USDS, tol),psi);
    end
    function [n] = nj_single(taui,nj,zetaj, psi, psimax, M, USDS, tol)
        %Helping function to nj_static
        %This function is needed becase integral() passes a phi as a vector
        %argument, while in turn integral() itselv can't handel vectors as
        %endpoints.
        % the varible lim is either +1 or -1.
        import Shock_pkg_new.Shock
        
        v0=real(sqrt(2*(psimax-psi)));
        n=integral(@(v) Shock.fj_static(taui,nj,zetaj, psi, M, v), -Inf,USDS*v0, 'RelTol', tol);
    end

end %end methods (static)

end %end class



