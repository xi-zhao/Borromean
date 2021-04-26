function varargout=dimerrf2(B,E,omega)
alp=1;
eta=40/6;
n1=500;
n2=500;
if imag(B)==0
    zeroi=0.01*complex(0,1);
else 
    zeroi=0;
end
[x1,wx1]=lgwt(n1,0,1);
[x2,wx2]=lgwt(n2,0,1);
k_perp=[x1;1./x2];
wk_perp=[wx1;wx2./x2.^2];
k_z=[-1./x2;-x1;x1;1./x2];
wk_z=[wx2./x2.^2;wx1;wx1;wx2./x2.^2];
[ndK_perp,ndK_z]=ndgrid(k_perp,k_z);
[WK_perp,WK_z]=ndgrid(wk_perp,wk_z);

if and(imag(B)==0,abs(B)<1) 
    Eth=real(-alp^2-B^2-2*alp*B);
elseif  and(abs(B)<alp,imag(B)~=0)
    Eth=real(-alp^2-B^2);
else
    Eth=real(-4*B);
end
x1=B./ndK_perp/alp;bkp1=-x1-sqrt(x1.^2+1);bkm1=-x1+sqrt(x1.^2+1);
R_kp_up=1./sqrt(1+bkm1.^2);
R_km_up=-1./sqrt(1+bkp1.^2);
R_kp_down=-bkm1./sqrt(1+bkm1.^2);
R_km_down=bkp1./sqrt(1+bkp1.^2);
L_kp_up=-R_km_down;
L_km_up=R_kp_down;
L_kp_down=-R_km_up;
L_km_down=-R_kp_up;



Kmod_square=ndK_perp.^2+ndK_z.^2;

psiRp=R_kp_up./(E+Eth-eta*Kmod_square-(Kmod_square-2*B+2*sqrt(alp^2*ndK_perp.^2+B^2)));
psiRm=R_km_up./(E+Eth-eta*Kmod_square-(Kmod_square-2*B-2*sqrt(alp^2*ndK_perp.^2+B^2)));
psiLp=conj(psiRp).*(-conj(R_km_down).*L_kp_up+conj(R_km_up).*L_kp_down)+conj(psiRm).*(conj(R_kp_down).*L_kp_up-conj(R_kp_up).*L_kp_down);
psiLm=conj(psiRp).*(-conj(R_km_down).*L_km_up+conj(R_km_up).*L_km_down)+conj(psiRm).*(conj(R_kp_down).*L_km_up-conj(R_kp_up).*L_km_down);

denomip=omega-(E+Eth-eta*Kmod_square-(Kmod_square-2*B+2*sqrt(alp^2*ndK_perp+B^2)))+zeroi;
denomim=omega-(E+Eth-eta*Kmod_square-(Kmod_square-2*B-2*sqrt(alp^2*ndK_perp+B^2)))+zeroi;

part1=-imag(sum(sum(ndK_perp.*WK_perp.*WK_z.*psiRp.*psiLp./denomip)));
part2=-imag(sum(sum(ndK_perp.*WK_perp.*WK_z.*psiRm.*psiLm./denomim)));
varargout{1}=part1+part2;
varargout{2}=part1;
varargout{3}=part2;
end


