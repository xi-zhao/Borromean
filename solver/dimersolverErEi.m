function res=dimersolverErEi(reE,imagE,kfasneg,eta,B,alp,cutoff)
% dimer bound state solver
E=real(reE)+complex(0,1).*real(imagE);

switch nargin
    case 4
        B=0;alp=1;invcutoff=0;
    case 5
        alp=1;invcutoff=0;
    case 6
        invcutoff=0;
    case 7
        invcutoff=1/cutoff;
end

na=200;

[x1,w1]=lgwt(na,0,1);
[x2,w2]=lgwt(na,invcutoff,1);
k_perp=[x1;1./x2];
wk_perp=[w1;w2./x2.^2];
k_z=k_perp;
wk_z=[w1;w2./x2.^2];

K_perp=repmat(k_perp,1,2*na);
K_z=repmat(k_z',2*na,1);

if and(imag(B)==0,abs(B)<1) 
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=real(-4*B);
end



Kmod_square=K_perp.^2+K_z.^2;

part1=pi*kfasneg/(eta+1);
x1=B./k_perp./alp;
bkp1=x1+sqrt(x1.^2+1);bkm1=x1-sqrt(x1.^2+1);
gamma_kp_up=bkp1./sqrt(1+abs(bkp1).^2);
gamma_km_up=bkm1./sqrt(1+abs(bkm1).^2);
gamma_kp_down=1./sqrt(1+abs(bkp1).^2);
gamma_km_down=1./sqrt(1+abs(bkm1).^2);


detU=gamma_kp_up.*gamma_km_down-gamma_kp_down.*gamma_km_up;
L_kp_up=gamma_kp_up;
L_km_up=gamma_km_up;
L_kp_down=gamma_kp_down;
L_km_down=gamma_km_down;

R_kp_up=gamma_km_down./detU;
R_km_up=-gamma_kp_down./detU;
R_kp_down=-gamma_km_up./detU;
R_km_down=gamma_kp_up./detU;




part2temp_1=2*K_perp/(eta+1)./Kmod_square;


part2temp_2=2*R_kp_up.*L_kp_up.*K_perp./(E+Eth-(eta+1)*Kmod_square-2*alp*sqrt(K_perp.^2+B^2)+2*alp*B)+...
    2*R_km_up.*L_km_up.*K_perp./(E+Eth-(eta+1)*Kmod_square+2*alp*sqrt(K_perp.^2+B^2)+2*alp*B);
part2temp=part2temp_1+part2temp_2;

part2=wk_perp'*part2temp*wk_z;
res=part1-part2;
end

