function cost=hyperspectral_cost(x,Y,AdditionalData)

MIE_A=x(1);
MIE_B=x(2);
Hbo=x(3);
hbr=x(4);
cco_diff=x(5);
WF=x(6);
oxredCCO=x(7);
oxCCO=x(8);
redCCO=x(9);
Fat=x(10);
DC=x(11);

ext_coef=AdditionalData.ext_coef;
lambda = AdditionalData.Lambda;
Reff = AdditionalData.Reff;
rho = AdditionalData.rho;


mua = log(10)*(1E-6*ext_coef(:,1)*Hbo + 1E-6*ext_coef(:,2)*hbr +1E-6*ext_coef(:,3)*cco_diff +...
    ext_coef(:,4)*WF*100 + 1E-6*ext_coef(:,5)*oxredCCO+1E-6*ext_coef(:,6)*oxCCO+ext_coef(:,7)*1E-6*redCCO+ext_coef(:,8)*Fat);
if(any(mua<0))
    cost=1E9;
    return
end

log(Y)=log(A)-B*log(lamda)

mus = MIE_A*lambda.^(-MIE_B);
mue = (3.*mua.*(mua+mus)).^(1/2)-DC;
y=mue;
D   = 1./(3*mua+3*mus);
z0  = 1./(mua+mus);
zb  = ((1+Reff)/(1-Reff))*2*D;
x0 = z0.^2 + rho^2;
x1 = (z0+2*zb).^2 + rho^2; 

y = -log((1./4*pi).*((z0.*(mue+1./sqrt(x0)).*exp(-mue.*sqrt(x0))./x0) +...
    ((z0 + 2*zb).*(mue + 1./sqrt(x1)).*exp(-mue.*sqrt(x1))./x1)))-DC;

cost=mean((Y-y).^2);