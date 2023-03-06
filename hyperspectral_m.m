function [dx,y,varargout]=hyperspectral(t,x,u,...
    MIE_A,MIE_B,WF,oxredCCO,oxCCO,redCCO,Fat,DC,...
    gain_hbo,gain_hb,gain_diffcco,varargin)

varargin=varargin{1};
% x = [HbO, Hb, CCO_redux]
lambda = varargin{1}.Lambda;
ext_coef = varargin{1}.ext_coef;
rho =varargin{1}.rho;
Reff = varargin{1}.Reff;

%Reff =0.493;
% rho
% WF
% CCO
% MIE_A
% MIE_B
% 
% gain_hbo
% gain_hb
% gain_diffcco




mua = log(10)*(1E-6*ext_coef(:,1)*x(1) + 1E-6*ext_coef(:,2)*x(2) +1E-6*ext_coef(:,3)*x(3) +...
    ext_coef(:,4)*WF + 1E-6*ext_coef(:,5)*oxredCCO+1E-6*ext_coef(:,6)*oxCCO+ext_coef(:,7)*1E-6*redCCO+ext_coef(:,8)*Fat+DC);
mus = MIE_A*lambda.^(-MIE_B/1000);
mue = (3.*mua.*(mua+mus)).^(1/2);
D   = 1./(3*mua+3*mus);
z0  = 1./(mua+mus);
zb  = ((1+Reff)/(1-Reff))*2*D;
x0 = z0.^2 + rho^2;
x1 = (z0+2*zb).^2 + rho^2; 

y = (1./4*pi).*((z0.*(mue+1./sqrt(x0)).*exp(-mue.*sqrt(x0))./x0) +...
    ((z0 + 2*zb).*(mue + 1./sqrt(x1)).*exp(-mue.*sqrt(x1))./x1));



dx = [gain_hbo*u; gain_hb*u; gain_diffcco*u];

if(nargout>2)
    % Return dF/dx
    zeros(3);
end

if(nargout>3)
    % Return dy/dx
    %-->TODO
end


