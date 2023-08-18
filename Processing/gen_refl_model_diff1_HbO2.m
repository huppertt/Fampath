function [diff_1_model_refl] = gen_refl_model_diff1(parameters,xdata,WF,HHb)
%GEN_REFL_MODEL This function generates the reflectance model for
%wavelenghts from 680-825nm 
%   This function takes inputs and uses them to generate a reflectance
%   model for wavelengths from 680-850nm. The equation is taken from Kienle
%   et al. 2 light intensities, z0 and zb are needed to fulfil the
%   extrapolated boundary conditions. 
coeff = load_extinction_coeff;

rho = 30; % source-detector separation 

% parameters = [WF, HHb, HbO2, a, b] 

WF = parameters(1);
HHb = parameters(2);
HbO2 = parameters(3);
a = parameters(4);
b = parameters(5);

indices = (xdata-650+1); %indices to read off the table
mua = log(10).*HbO2.*coeff(indices,3) + log(10).*HHb.*coeff(indices,2) + WF.*coeff(indices,4); 
mus = a.*xdata'.^(-b); % reduced scattering coefficient from Matcher et al.
mue = (3.*mua.*(mua + mus)).^1/2; % effective attenuation coefficient
D = 1./(3.*(mua+mus));

z0 = 1./(mua+mus);
zb = ((1+0.493)/(1-0.493))*2.*D; %0.493 is Reff at n=1.4

x0 = z0.^2 + rho^2;
x1 = (z0+2.*zb).^2 + rho^2; 

model_refl = (1./4*pi).*((z0.*(mue+1./sqrt(x0)).*exp(-mue.*sqrt(x0))./x0) + ((z0 + 2*zb).*(mue + 1./sqrt(x1)).*exp(-mue.*sqrt(x1))./x1));
model_refl = smooth(model_refl); 

diff_1_model_refl = diff(smooth(model_refl,1)); 
diff_2_model_refl = diff(smooth(diff_1_model_refl)); 
end


