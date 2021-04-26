[xi,wi]=lgwt(30,-1,1);
funold=cos(xi).*sin(xi);
x=0:0.05:2;
dimx=length(x);
%generating x-dim
ximxj=repmat(xi',[length(xi),1])-repmat(xi,[1,length(xi)]);
ximxj=ximxj+eye(length(xi));
ximxjinv=prod(1./ximxj);
xmxj=prod(repmat(x',[1,length(xi)])-repmat(xi',[length(x),1]),2);
fun=((xmxj*ximxjinv)./(repmat(x',[1,length(xi)])-repmat(xi',[dimx,1])))*funold;
fun=sum(fun,2);
plot(x,fun)
hold on 
plot(xi,funold,'.')
