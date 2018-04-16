function [w_out] = w_approx(eps,Upsilon)

%if nargin <= 2
%    threshold=1e-1;
%end

%eps=reshape(eps,1,[]);
w_out=zeros(size(eps));

% Look-up table for the value of w at the positive threshold
% end                 w(end)
% 1.000000000000000   1.475633535051814

threshold=0.1;
Ip1=find(eps>0 & eps<=threshold);
%Ip1=find(eps>0);
Ep1=eps(Ip1);
w_out(Ip1)=wp(Ep1);

Ip2=find(eps>threshold);
Ep2=eps(Ip2);
wp_asymp=@(E) (E- (1-E).^3/48 )*sqrt(2);
%                                    vvv -- value of w at E=1
w_out(Ip2)=wp_asymp(Ep2)-wp_asymp(1)+1.475633535051814;
%w_out(Ip2)=wp_asymp(Ep2)-wp_asymp(threshold)+0.272464896153391;
%w_out(Ip2)=(Ep2-min(Ep2))*sqrt(2) + wp(max(Ep1));


%  Look-up table for the value of w at the first negative threshold
%  threshold           w(threshold)
%  -0.2969   -0.4146
% -0.961207596286194  -1.163402271526147
threshold=-0.961207596286194;
In1=find(eps<0 & eps>=threshold);
En1=eps(In1);
w_out(In1)=wn(En1);

In2=find(eps<-0.5);
En2=eps(In2);
% Intermediary approximation of w.
% The numerical factor is 1/sqrt(Shock_pkg.Feps(-1.8))
%wn_approx=@(E) E*1.0427;%/sqrt(Shock_pkg.Feps(-.4));
wn_approx=@(E) E*1.034099637139201; %/sqrt(Shock_pkg.Feps(threshold-1))
w_out(In2)=wn_approx(En2)-wn_approx(threshold)-1.163402271526147;
%wn_asym=@(E) E;
%w_out(In2)=wn_asym(En2)-wn_asym(max(En2)) + wn(min(En1));



threshold=-2;
In3=find(eps<threshold);
En3=eps(In3);
wn_asymp=@(E) E; %-1./E;
%                                      vvv -- value of w at E=-50
w_out(In3)=wn_asymp(En3)-wn_asymp(-50)-50.330999102857149;

Iz=find(eps==0);
%Ez=eps(Iz);
w_out(Iz)=zeros(size(Iz));

w_out=Upsilon*w_out;



end

function Wp = wp(Ep)
sqrl=sqrt(log(8./Ep));
Wp=Ep.*sqrl/sqrt(2)+2*sqrt(2*pi)*erfc(sqrl);
end

function Wn = wn(En)
sqrl=sqrt(log(-8./En));
Wn=En.*sqrl/sqrt(2)-2*sqrt(2*pi)*erfc(sqrl);
end

