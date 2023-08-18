function out = sim_refl_model(parameters,wavelengths,options,returntype)
%   This function takes inputs and uses them to generate a reflectance
%   model for wavelengths from 680-850nm. The equation is taken from Kienle
%   et al. 2 light intensities, z0 and zb are needed to fulfil the
%   extrapolated boundary conditions. 

if(nargin<4)
    returntype='all';
end

a = parameters(1);
b = parameters(2);
scale=parameters(3);

mua = log(10)*parameters(:,4:end)*options.coeff';

mus = a.*(wavelengths/500).^(-b); % reduced scattering coefficient from Matcher et al.
mue = (3.*mua.*(mua + mus)).^1/2; % effective attenuation coefficient
D = 1./(3.*(mua+mus));

z0 = 1./(mua+mus);
zb = ((1+options.Reff)/(1-options.Reff))*2.*D; %0.493 is Reff at n=1.4

x0 = z0.^2 + options.rho^2;
x1 = (z0+2.*zb).^2 + options.rho^2; 

model.reflectance =(1./4*pi).*((z0.*(mue+1./sqrt(x0)).*exp(-mue.*sqrt(x0))./x0) +...
    ((z0 + 2*zb).*(mue + 1./sqrt(x1)).*exp(-mue.*sqrt(x1))./x1));

model.reflectance=scale*exp(-(mua+mus));

model.reflectance = options.smoother(model.reflectance); 
model.diff_1_reflectance = options.smoother(diff(model.reflectance,1)); 
model.diff_2_reflectance = options.smoother(diff(model.diff_1_reflectance)); 

if(strcmp(lower(returntype),'all'))
     out=[model.reflectance,...
     model.diff_1_reflectance,...
     model.diff_2_reflectance];
elseif(strcmp(lower(returntype),'refl'))
    out=model.reflectance;
elseif(strcmp(lower(returntype),'diff_1'))
    out=model.diff_1_reflectance;
elseif(strcmp(lower(returntype),'diff_2'))
    out=model.diff_2_reflectance;
else
    out=[];
end