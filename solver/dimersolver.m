function res=dimersolver(E,kfasneg,eta,B,alp,cutoff)
% dimer bound state solver

switch nargin
    case 3
        B=0;alp=1;invcutoff=0;
    case 4
        alp=1;invcutoff=0;
    case 5
        invcutoff=0;
    case 6
        invcutoff=1/cutoff;
end

na=20;

[x1,w1]=lgwt(na,0,1);
[x2,w2]=lgwt(na,invcutoff,1);
k_perp=[x1;1./x2];
wk_perp=[w1;w2./x2.^2];
k_z=[-1./x2;-x1;x1;1./x2];
wk_z=[w2./x2.^2;w1;w1;w2./x2.^2];

[K_perp,K_z]=ndgrid(k_perp,k_z);
[WK_perp,WK_z]=ndgrid(wk_perp,wk_z);


if and(imag(B)==0,abs(B)<1) 
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=real(-4*B);
end



Kmod_square=K_perp.^2+K_z.^2;

part1=pi*kfasneg/(eta+1);

xi=@(k_perp,k_z,sigma) k_perp.^2+k_z.^2-2*B+2*sigma*sqrt(alp^2*k_perp.^2+B^2);
epsilonb=@(k_perp,k_z) eta*(k_perp.^2+k_z.^2);

b=@(k,sigma) (B+sigma*sqrt(B^2+alp^2*k.^2))./(alp*k);
gama=@(k,sigma,arrow) (arrow==1)*b(k,sigma)./sqrt(1+abs(b(k,sigma)).^2)+(arrow==-1)./sqrt(1+abs(b(k,sigma)).^2);
L=@(k,sigma) gama(k,sigma,1);
R=@(k,sigma) sigma*gama(k,-sigma,-1)./(gama(k,1,1).*gama(k,-1,-1)-gama(k,-1,1).*gama(k,1,-1));

O=@(sigma) K_perp.*(1./(2*(1+eta)*Kmod_square)+R(K_perp,sigma).*L(K_perp,sigma)./(E+Eth-xi(K_perp,K_z,sigma)-epsilonb(K_perp,K_z)));
part2= sum((O(1)+O(-1)).*WK_perp.*WK_z,'All');
res=part1-part2;


end

