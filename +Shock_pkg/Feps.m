function [F] = Feps(e)

k=1-e;


Ip=find(e>0);
kp=k(Ip).^2;
[Kp,Ep]=ellipke(kp);
Fp=(Ep./Kp+kp-1)./kp;

Ulim=0.999;
IU=find(e>Ulim);
FU=.5-k(IU).^2/16;


In=find(e<=0);
kn=k(In).^(-2);
[Kn,En]=ellipke(kn);
Fn=En./Kn;

F=ones(size(e));
F(Ip)=Fp;
F(In)=Fn;
F(IU)=FU;

%{
if e>0
    ks=k.^2;
    [K,E]=ellipke(ks);
    f=(E./K+ks-1)./ks;
else
    ks=k.^(-2);
    [K,E]=ellipke(ks);
    f=E./K;
end
%}
end

