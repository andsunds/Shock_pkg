classdef Shock_MB < Shock_pkg.Shock
% Shocks using the assumption that the electrons are Maxwell-Boltzmann
% distributed.
%
% Author: Andréas Sundström (c)
%

properties
    % This class has no own properties other than what is inhereted from
    % the Shock superclass.
end


methods
    function obj = Shock_MB(m,Z,n, tau, speed_type,speed, phimaxmin_in, tol)
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
        % If everything is ok, this finds and sets phimax
        if ( obj.ion_species>0 )&&( ~isempty(obj.V) )
            obj.phimax=find_phimax(obj, phimaxmin_in(1));
            % Checks that the solution is physical (actual maximum at
            % phi=phimax).
            if obj.charge_dens(+1, obj.phimax)<0 || obj.charge_dens(-1, obj.phimax)<0
                fprintf('ERROR: found phimax = %1.4f is NOT a true maximum.\n',obj.phimax)
                obj.phimax=NaN;
                return
            elseif length(phimaxmin_in)==2
                phimin_in=phimaxmin_in(2);
                obj.phimin=find_phimin(obj, phimin_in);
            else
                phimin_in=obj.phimax * 0.5;
                obj.phimin=find_phimin(obj, phimin_in);
            end
            if obj.phimin<0
                %obj.phimax=NaN;
                obj.phimin=NaN;
            end
        end
    end %end Constructor
    
    function [n_el] = ne(obj, phi)
        % The total elctron density, due to ALL ion species.
        n_el=0; %init
        for i=1:obj.ion_species %summing over all species
        n_el=n_el+obj.ne_static(obj.m(i),obj.Z(i),obj.n(i),...
            phi, obj.phimax, obj.V, obj.tau);
        end
    end

end %end methods

methods (Access=protected)
    function [phimax]=find_phimax(obj, phimax_in)
        % Finds the correct value of phimax, used in the constructor.
        phimaxin=phimax_in;%*(obj.Mach^2*obj.tau/2); %initial guess
        %The Sagdeev potential has to be 0 at phi=phimax
        phimax=fzero(@(phim) obj.Phi_static(+1,obj.m,obj.Z,obj.n,...
                        phim, phim, obj.V, obj.tau, obj.tol), phimaxin,...
                        optimset('tolfun',obj.tol));
    end
    
    function [phimin]=find_phimin(obj, phimin_in)
        % Finds phimin is obj is a valid shock.
        %Checks to see if all is well
        if isnan(obj.phimax)
            phimin=NaN;
        else
            %The DS Sagdeev potential has to be 0 at phi=phimmin
            phimin=fzero(@(phim) obj.Phi_static(-1,obj.m,obj.Z,obj.n,...
                            phim, obj.phimax, obj.V, obj.tau, obj.tol), phimin_in,...
                            optimset('tolfun',obj.tol));
            if phimin>=obj.phimax
                fprintf('The specified phimin_in found the wrong root of Phi.\nNow trying in the interval [0, 1-tol]*phimax.\n\n')
                %Now we instead define an interval in which to look
                phimin_in=[0,(1-obj.tol)]*obj.phimax;
                if obj.Phi(-1, phimin_in(2) )*obj.Phi(-1, 0)<0
                    phimin=find_phimin(obj, phimin_in);
                else
                    phimin=NaN;
                    fprintf('ERROR: No phimin in the interval [0, 1-tol]*phimax.\n')
                end
            end
        end
    end
 
end %end methods (proteccted)

methods (Static=true, Access=protected)
    % Defining static functions which in turn are defined to be used in the
    % above functions, these should and can not be used out side this class
        
    function [ne_j] = ne_static(mj,Zj,nj, phi, phimax, V, tau)
        % The elctron density, due to ONE ion species. To get the total
        % value, sum this function over all ions. This is assuming that the
        % electron distribution function is walways a pure
        % Maxwell-Boltzmann. 
        % The ion density (of this specific species) far upstream from the
        % shock, but still containing reflected ions. ("n" is without the
        % refelcted ions.) 
        njInfUS=0.5*nj*( 1 + 2*erf(V*sqrt(mj/2)) +...
            erf(real(sqrt(2*Zj*phimax/mj)-V)*sqrt(mj/2)) );
        
        % Given a pure M-B, the distribution integrates up to this density:
        ne_j=Zj*njInfUS*exp(phi/tau);
    end
    
    function [rho] = charge_dens_static(USDS, m,Z,n, phi, phimax, V, tau, tol)
        %  Note reqiures that m, Z, and n are all the same length.
        import Shock_pkg.Shock_MB
        L=length(Z);
        rho=0;
        for j=1:L
        rho=rho+( -Shock_MB.ne_static(m(j),Z(j),n(j), phi, phimax, V, tau)...
                  +Z(j)*Shock_MB.nj_static(USDS,m(j),Z(j),n(j), phi, phimax, V, tol) );
        end
    end

    function [Phi_single] = Phi_static(USDS,m,Z,n, phi, phimax, V, tau, tol)
        %Must be static to be able to use this in find_phimax().
        import Shock_pkg.Shock_MB
        if USDS==1
            Phi_single=integral(@(phiP) Shock_MB.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, V, tau, tol), 0, phi, 'RelTol',tol);
        elseif USDS==-1
            Phi_single=integral(@(phiP) Shock_MB.charge_dens_static(...
                USDS,m,Z,n, phiP, phimax, V, tau, tol), phimax, phi, 'RelTol',tol);
        else
            Phi_single=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
        end
    end
end %end methods



end %end class

