function [ Fp,Fn ] = Flux( k,phi,r )
% Fp(k,phi,r) = integral(Pr(phi-t)*exp(ikt),-pi/2,pi/2);
k1 = k; k = abs(k); % deal with k<0
Fp = r^(-k)*exp(1i*k*phi)/(2*pi*1i)*(log(1+1i*r*exp(-1i*phi))-log(1-1i*r*exp(-1i*phi)))...
    -r^k*exp(1i*k*phi)/(2*pi*1i)*(log(1-1i*r*exp(1i*phi))-log(1+1i*r*exp(1i*phi)))...
    +r^k/2*exp(1i*k*phi);
for m = 0:(k-1)/2
   Fp = Fp+exp(1i*(k-1)*phi)/pi*(r^(k-2*m-1)-r^(2*m+1-k))/(2*m+1)*exp(1i*m*(pi-2*phi)); 
end
% Fn(k,phi,r) = integral(Pr(phi-t)*exp(ikt),Re(exp(it))<0)
Fn = r^k*exp(1i*k*phi)-Fp;
if(k1<0)
    Fp = conj(Fp); Fn = conj(Fn);
elseif(k1==0)
    Fp = real(Fp); Fn = real(Fn);
end

end

