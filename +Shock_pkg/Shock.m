classdef (Abstract=true) Shock < matlab.mixin.CustomDisplay
% Base class for Shocks. This is an abstrac class.
% This class contains some of the basic properties and methods needed to
% represent the 1D electrostatic shocks in plasmas studied in this project.
%
% The inheritance from matlab.mixin.CustomDisplay, is just so that the
% property display can be changed.
%
% Author: Andréas Sundström (c)
%
    
properties (SetAccess = protected)
    tol  %The tolerance to which the integrations should be carried
    
    %These properties can be vectors, but they have to be the same length.
    n    % ion densities (normalized so that main ion density is 1)
    Z    % charges
    m    % masses
    %T   % ion temperatures. Not implemented!
    %No longer vector quantities
    ion_species % the number of ion species
    tau         % electron-ion temp ration
    V           % flow speed
    phimax      % Max value of the electrostatic potential
    phimin % Min value of the electrostatic potential in downstream region
    
end

properties (Dependent)
    Mach   %Mach number
    cs     %Speed of sound
    F      %phimax=F*(Mach^2*tau/2)
    
    lambda %Wavelength of DS oscillations
end

    
methods
    function obj = Shock(m,Z,n, tau, speed_type, speed, tol)
        %Constructor
        if nargin==0 %this is due to the way you call a superclass constructor
            %Do nothing
        else
            %Set values upon construction.
            obj.tol=tol;
            obj.tau=tau;   %electron-ion temp ration
            obj.n=n/n(1);  %ion densities (normalized so that main ion density is 1)
            obj.Z=Z;       %charges
            obj.m=m;       %masses
            
            %The variable speed_type determines which type of measurement
            %speed is (Mach or V).
            if strcmp(speed_type,'Mach')
                obj.V=speed*obj.cs;
            elseif speed_type == 'V'
                obj.V=speed;
            else
                obj.V=[]; %obj.Mach=[];
                fprintf('ERROR: the speed type must be either "Mach" och "V".\n')
            end
        end
        %checking so that m, n, and Z are the same length
        L1=length(obj.m);L2=length(obj.Z);L3=length(obj.n);
        if (L1==L2)&&(L2==L3)
            obj.ion_species=L1;
        else
            obj.n=[]; obj.Z=[]; obj.m=[];
            fprintf('ERROR: m, Z, and n are not the same length!\n\n')
            return
        end
    end %end consructor
    
    function cs=get.cs(obj)
        % Get function for the speed of sound (dependent variable)
        cs=sqrt(obj.Z(1)*obj.tau/obj.m(1));
    end
    function M=get.Mach(obj)
        % Get function for the Mach number(dependent variable)
        %cs=sqrt(obj.Z(1)*obj.tau/obj.m(1));
        M=obj.V/obj.cs;
    end
    function Fval = get.F(obj)
        %Get function of for the F vaule
        Fval=[obj.phimax,obj.phimin]/(obj.Mach^2*obj.tau/2);
    end
    function L=get.lambda(obj)
        %Get function for lambda
        L = sqrt(2)*integral(@(p)1./sqrt(-obj.Phi(-1,p)),...
            obj.phimin, obj.phimax, 'RelTol', obj.tol);
    end

    

    function [X,phi,E,rho] = find_phi(obj,Xmin,Xmax)
        % This function calcualtes the potential phi as a function of x,
        % also the electric field E=dphi/dx comes "free of charge".
        
        % Functions specifying the derivatives of phi (G)
        odefunUS=@(x, G) [G(2);-obj.charge_dens(+1, G(1))];
        odefunDS=@(x, G) [G(2);-obj.charge_dens(-1, G(1))];
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
            X  = [flip(xDS,1);xUS];
            phi= [flip(phiDS(:,1),1);phiUS(:,1)];
            E  =-[flip(phiDS(:,2),1);phiUS(:,2)];
            
            if nargout>=4
                rho=zeros(length(X),1);
                for i=1:length(rho)
                    if X(i)>0
                        dG=odefunUS(X(i), [phi(i),E(i)]);
                    else
                        dG=odefunDS(X(i), [phi(i),E(i)]);
                    end
                    rho(i)=-dG(2);
                end
            end
        else
            X=[];phi=[];E=[];rho=[];
            fprintf('ERROR: not a maximum at phi=phimax, charge density is negative.\n')
        end
    end

    function [ni] = nj(obj, USDS, phi, index)
        %Calculates the UpStream ion density of ion species number 'index'.
        %This is done using the static function nj_single.
        % The argument USDS must be either +1 (US) or -1 (DS).
        if abs(USDS)==1
            ni=arrayfun(@(x) obj.nj_single(obj.m(index),obj.Z(index),...
                obj.n(index), x, obj.phimax, obj.V, USDS, obj.tol), phi);
        else
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            ni=NaN;
        end
    end
    

    function [rho] = charge_dens(obj, USDS, phi)
        % Retruns the charge density either US or DS at potential phi .
        % The argument USDS must be either +1 (US) or -1 (DS).
        if abs(USDS)==1
            rho=0;
            for i=1:obj.ion_species
                rho=rho+obj.Z(i)*obj.nj(USDS, phi, i);
            end
            rho=rho-obj.ne(phi);
        else 
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            rho=NaN;
        end
    end

    function [PHI] = Phi(obj, USDS, phi)
        % Calculates the Sagdeev potential.
        % As per the theory, the Sagdeev potential is just the chagre
        % density integrated over the electrostatic potential. This is used
        % to calculate the maximum value of phi. 
        %    The argument USDS must be either +1 (US) or -1 (DS).
        if abs(USDS)==1 
            PHI=arrayfun(@(x) obj.Phi_single(USDS,x),phi);
        else
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
            PHI=NaN;
        end
    end
end  %end methods

methods (Access=protected)
    function propgrp = getPropertyGroups(~)
        %Function defining how to display the properties
        proplist = {'tol','n','Z','m','ion_species','tau','Mach','V','cs',...
            'phimax','phimin','F'};
        propgrp = matlab.mixin.util.PropertyGroup(proplist);
    end
    function [PHI] = Phi_single(obj, USDS, phi)
        %Function that takes a single phi value, used in Phi().
        if USDS==1
            PHI=integral(@(phiP) obj.charge_dens(USDS, phiP), 0, phi, 'RelTol',obj.tol);
        elseif USDS==-1
            PHI=integral(@(phiP) obj.charge_dens(USDS, phiP), obj.phimax, phi, 'RelTol',obj.tol);
        else
            PHI=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
        end
    end
end


methods (Abstract=true, Static=true, Access=protected)
    [ne_j] = ne_static(mj,Zj,nj, phi, phimax, V, tau);
    [rho] = charge_dens_static(USDS, m,Z,n, phi, phimax, V, tau, tol);
    [Phi_single] = Phi_static(USDS,m,Z,n, phi, phimax, V, tau, tol);
end



methods (Static=true, Access=protected)
    % Defining static functions which in turn are used in the definitions
    % of some of the above functions, these should and can not be used
    % outside of this class or its subclasses. 
    
    function [dist] = fj_static(mj,Zj,nj, phi, V, v)
        % The ion distributin function of a specific ion species
        dist = nj*sqrt(mj/(2*pi)) * exp(-(mj/2)*(sqrt(v.^2+2*Zj*phi/mj)-V).^2);
    end
    function [n] = nj_static(USDS,mj,Zj,nj, phi, phimax, V, tol)
        %Calclulats the upstream or downstream ion density for ONE ion
        %species, due to the value of the potential and the maximum value
        %of the potential. 
        import Shock_pkg.Shock
        n=arrayfun(@(x) Shock.nj_single(mj,Zj,nj, x, phimax, V, USDS, tol), phi);
    end
    function [n] = nj_single(mj,Zj,nj, phi, phimax, V, USDS, tol)
        %Helping function to nj_static
        %This function is needed becase integral() passes a phi as a vector
        %argument, while in turn integral() itselv can't handel vectors as
        %endpoints.
        % the varible lim is either +1 or -1.
        import Shock_pkg.Shock
        v0=real(sqrt(2*Zj*(phimax-phi)/mj));
        n=integral(@(v) Shock.fj_static(mj,Zj,nj, phi, V, v), -Inf,USDS*v0, 'RelTol', tol);
    end

end %end methods (static)

end %end class



