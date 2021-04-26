function res=dimersolvercoeffm(E,kfasneg,eta,B,alp,cutoff)

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

m1=10;
m2=10;
[x1,w1]=lgwt(m1,0,1);
[x2,w2]=lgwt(m2,invcutoff,1);
[theta,~]=lgwt(60,0,2*pi);
k_perp=[x1;1./x2];
wk_perp=[w1;w2./x2.^2];
k_z=[-1./x2;-x1;x1;1./x2];
wk_z=[w2./x2.^2;w1;w1;w2./x2.^2];
[K_perpd,K_zd]=ndgrid(k_perp,k_z);
[Wk_perpdd,Wk_zd]=ndgrid(wk_perp,wk_z);
[K_perp,K_z,Q_perp,Q_z]=ndgrid(k_perp,k_z,k_perp,k_z);
[~,~,WQ_perp,WQ_z]=ndgrid(wk_perp,wk_z,wk_perp,wk_z);
if and(imag(B)==0,abs(B)<1)
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=real(-4*B);
end
E=E+Eth;
D=pi*kfasneg/(1+eta)-sum(K_perpd./(K_perpd.^2+K_zd.^2)/(1+eta).*Wk_perpdd.*Wk_zd,'all');
D=D*eye(2*length(k_perp)*length(k_z));

xi=@(k_perp,k_z,sigma) k_perp.^2+k_z.^2-2*B+2*sigma.*sqrt(alp^2*k_perp.^2+B^2);
epsilonb=@(k_perp,k_z) eta*(k_perp.^2+k_z.^2);

b=@(k,sigma) (B+sigma*sqrt(B^2+alp^2*k.^2))./(alp*k);
gama=@(k,sigma,arrow) (arrow==1)*b(k,sigma)./sqrt(1+abs(b(k,sigma)).^2)+(arrow==-1)./sqrt(1+abs(b(k,sigma)).^2); 
L=@(k,sigma) gama(k,sigma,1);
R=@(k,sigma) sigma*gama(k,-sigma,-1)./(gama(k,1,1).*gama(k,-1,-1)-gama(k,-1,1).*gama(k,1,-1));

numerator=@(sigmap,sigma) 1i*R(K_perp,sigma).*L(Q_perp,sigmap);
denominator=@(sigmap,sigma) E-xi(K_perp,K_z,sigma)-epsilonb(K_perp,K_z);


OffDiag=@(sigmap,sigma) Q_perp./(pi).*numerator(sigmap,sigma)./denominator(sigmap,sigma);
Re=@(x) reshape(x.*WQ_perp.*WQ_z,length(k_perp)*length(k_z),[]);
O=[Re(OffDiag(1,1)),Re(OffDiag(1,-1));Re(OffDiag(-1,1)),Re(OffDiag(-1,-1))];

H=D+O;
res=eigs(H,1,'sm');


end