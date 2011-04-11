function [sig,sig2,osig,msig,tsig,tau]=avar(y,tau0)
% function [sig,sig2,osig,msig,tsig,tau]=avar(y,tau0)
% INPUTS:
%  y = signal
%  tau0 = sampling period (s)
%
% OUTPUTS:
% sig =  N samples STD DEV
% sig2 = Normal Allan STD DEV, 2 samples STD DEV.
% osig = Sigma(y)(tau) = Allan Standard Deviation with Overlapping estimate
% msig = Modified Allan Standard Deviation
% tsig = Time Allan Standard Deviation
% tau = measurement time (s).
%
% Copyright Alaa MAKDISSI 2003
% free for personal use only.
% Obtained from 
% http://www.alamath.com/index.php?option=com_content&task=view&id=19&Itemid=9

s=[];
x=[];
n=length(y);
jj=floor( log((n-1)/3)/log(2) );

for j=0:jj
    fprintf('.');
    m=2^j;
    tau(j+1)=m*tau0;
    D=zeros(1,n-m+1);
    for i=1:n-m+1
        D(i)=sum(y(i:i+m-1))/m;
    end
    %N sample
    sig(j+1)=std(D(1:m:n-m+1));
    %AVAR
    sig2(j+1)=sqrt(0.5*mean((diff(D(1:m:n-m+1)).^2)));
    %OVERAVAR
    z1=D(m+1:n+1-m);
    z2=D(1:n+1-2*m);
    u=sum((z1-z2).^2);
    osig(j+1)=sqrt(u/(n+1-2*m)/2);
   
    %MVAR
    u=zeros(1,n+2-3*m);
    for L=0:n+1-3*m
        z1=D(1+L:m+L);
        z2=D(1+m+L:2*m+L);
        u(L+1)=(sum(z2-z1))^2;
    end
   
    uu=mean(u);
    msig(j+1)=sqrt(uu/2)/m;
   
    %TVAR
    tsig(j+1)=tau(j+1)*msig(j+1)/sqrt(3);
   
end
