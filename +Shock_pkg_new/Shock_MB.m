classdef Shock_MB < Shock_pkg_new.Shock
% Shocks using the assumption that the electrons are Maxwell-Boltzmann
% distributed.
%
% Author: Andréas Sundström
%

properties
    % This class has no own properties other than what is inhereted from
    % the Shock superclass.
end


methods
    function obj = Shock_MB(Z,n, taui, Mach, psimaxmin_in, tol)
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
        
            obj.psimax=find_psimax(obj, psimaxmin_in(1));
            % Checks that the solution is physical (actual maximum at
            % phi=phimax).
            if obj.charge_dens(+1, obj.psimax)<0 || obj.charge_dens(-1, obj.psimax)<0
                fprintf('ERROR: found phimax = %1.4f is NOT a true maximum.\n',obj.psimax)
                obj.psimax=NaN;
                return
            elseif length(psimaxmin_in)==2
                psimin_in=psimaxmin_in(2);
                obj.psimin=find_psimin(obj, psimin_in);
            else
                psimin_in=obj.phimax * 0.5;
                obj.psimin=find_psimin(obj, psimin_in);
            end
            if obj.psimin<0
                %obj.phimax=NaN;
                obj.psimin=NaN;
            end
    end %end Constructor
    
    function [n_el] = ne(obj, psi)
        % The total elctron density, due to ALL ion species.
        n_el=0; %init
        for i=1 %summing over all species
        n_el=n_el+obj.ne_static(obj.taui, obj.Z(i),obj.n(i),...
            psi, obj.psimax, obj.M);
        end
    end

end %end methods

methods (Access=protected)
    function [phimax]=find_psimax(obj, psimax_in)
        % Finds the correct value of phimax, used in the constructor.
        psimaxin=psimax_in;%*(obj.Mach^2*obj.tau/2); %initial guess
        %The Sagdeev potential has to be 0 at phi=phimax
        phimax=fzero(@(psim) obj.Psi_static(+1,obj.Z,obj.n,...
                        psim, psim, obj.M, obj.taui, obj.tol), psimaxin,...
                        optimset('tolfun',obj.tol));
    end

    
    function [psimin]=find_psimin(obj, psimin_in)
        % Finds phimin is obj is a valid shock.
        %Checks to see if all is well
        if isnan(obj.psimax)
            psimin=NaN;
        else
            %The DS Sagdeev potential has to be 0 at phi=phimmin
            psimin=fzero(@(psim) obj.Psi_static(-1,obj.Z,obj.n,...
                            psim, obj.psimax, obj.M, obj.taui, obj.tol), psimin_in,...
                            optimset('tolfun',obj.tol));
            if psimin>=obj.psimax
                fprintf('The specified phimin_in found the wrong root of Phi.\nNow trying in the interval [0, 1-tol]*phimax.\n\n')
                %Now we instead define an interval in which to look
                psimin_in=[0,(1-obj.tol)]*obj.psimax;
                if obj.Psi(-1, psimin_in(2) )*obj.Psi(-1, 0)<0
                    psimin=find_psimin(obj, psimin_in);
                else
                    psimin=NaN;
                    fprintf('ERROR: No phimin in the interval [0, 1-tol]*phimax.\n')
                end
            end
        end
    end

end %end methods (proteccted)

methods (Static=true, Access=protected)
    % Defining static functions which in turn are defined to be used in the
    % above functions, these should and can not be used out side this class
        
    function [ne_j] = ne_static(taui, Zj,nj, psi, psimax, M)
        % The elctron density, due to ONE ion species. To get the total
        % value, sum this function over all ions. This is assuming that the
        % electron distribution function is walways a pure
        % Maxwell-Boltzmann. 
        % The ion density (of this specific species) far upstream from the
        % shock, but still containing reflected ions. ("n" is without the
        % refelcted ions.) 
        njInfUS=0.5*nj*( 1 + 2*erf(M*sqrt(taui/2)) +...
            erf(real(sqrt(2*psimax)-M)*sqrt(taui/2)) );
        % Given a pure M-B, the distribution integrates up to this density:
        ne_j=Zj*njInfUS*exp(psi);
    end
    
    function [rho] = charge_dens_static(USDS, Z,n, psi, psimax, M, taui, tol)
        %  Note reqiures that m, Z, and n are all the same length.
        import Shock_pkg_new.Shock_MB
        L=length(Z);
        rho=0;
        for j=1:L
        rho=rho+( -Shock_MB.ne_static(taui,Z(j),n(j), psi, psimax, M)...
                  +Z(j)*Shock_MB.nj_static(USDS,taui,n(j), psi, psimax, M, tol) );
        end
    end
   

    function [Psi_single] = Psi_static(USDS,Z,n, psi, psimax, M, taui, tol)
        %Must be static to be able to use this in find_phimax().
        import Shock_pkg_new.Shock_MB
        if USDS==1
            Psi_single=integral(@(psiP) Shock_MB.charge_dens_static(...
                USDS,Z,n, psiP, psimax, M, taui, tol), 0, psi, 'RelTol',tol);
        elseif USDS==-1
            Psi_single=integral(@(psiP) Shock_MB.charge_dens_static(...
                USDS,Z,n, psiP, psimax, M, taui, tol), psimax, psi, 'RelTol',tol);
        else
            Psi_single=[];
            fprintf('ERROR:the argument USDS must be either +1 (US) or -1 (DS).\n');
        end
    end
end %end methods



end %end class

