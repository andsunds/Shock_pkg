classdef Shock_Ad < Shock_pkg.Shock_tr
% Shocks using the assumption that the trapped electrons are adabatic and
% therefore have flat distribution function in the trapped regions.
%
% 
%
% Author: Andréas Sundström (c)
%


properties
    % This class has no own properties other than what is inhereted from
    % the Shock_tr superclass.
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
        proplist = {'tol','n','Z','m','tau','Mach','V','cs',...
            'phimax','phimin', 'trapping_coef'};
        propgrp = matlab.mixin.util.PropertyGroup(proplist);
    end
    function FVAL = zerofun_phiminmax(obj, phim)
        % The system of functions which we want to find the root of.
        FVAL=[obj.Phi_static(+1,obj.m,obj.Z,obj.n, phim(1), phim(1), phim(2), obj.V, obj.tau, obj.tol);
              obj.Phi_static(-1,obj.m,obj.Z,obj.n, phim(2), phim(1), phim(2), obj.V, obj.tau, obj.tol)];
    end

end %end methods


end %end class