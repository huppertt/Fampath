function hyperspectral_fit(filename)

data=load(filename);


% pathlength method
PL_method = 2;
% 1 = DPF, 2 = water fitting
%for water fitting
WaterFraction = 0.85;

%raw = Spectra; 

% Read in calibration file
%wl_cal = Wavelengths;
%wl = 780:900;
% Extinction coefficients
Ext=tissue_specific_extinction_coefficient_650to1000;

Ext_int=zeros(length(data.Wavelengths),size(Ext,2));
Ext_int(:,1)=data.Wavelengths;
for i=2:size(Ext,2)
    Ext_int(:,i)=interp1(Ext(:,1),Ext(:,i),Ext_int(:,1));
end

% HbO2, HHb, CCO
%ext_coeffs_inv = pinv(ext_coeffs);

DPF_dep=DPF_Lambda_Dependency_740to915;
DPF_dep_int=zeros(length(data.Wavelengths),size(DPF_dep,2));
DPF_dep_int(:,1)=data.Wavelengths;
for i=2:size(DPF_dep,2)
    DPF_dep_int(:,i)=interp1(DPF_dep(:,1),DPF_dep(:,i),DPF_dep_int(:,1),'nearest','extrap');
end

% Compute the reflectance data
reflect = data.Spectra./(100000*ones(size(data.Spectra,1),1)*data.Ref);
reflect_sm = sgolayfilt(reflect,4,45);
reflect_diff1 = sgolayfilt(diff(reflect_sm')',4,45);
reflect_diff2 = sgolayfilt(diff(reflect_diff1')',4,45);


mua = log10*HbO2*coeff(indices,3) + log(10).*HHb.*coeff(indices,2) + WF.*coeff(indices,4); 
mus = a.*lambda.^(-b); % reduced scattering coefficient from Matcher et al.
mue = (3.*mua.*(mua + mus)).^1/2; % effective attenuation coefficient
D = 1./(3.*(mua+mus));
z0 = 1./(mua+mus);
zb = ((1+0.493)/(1-0.493))*2.*D; %0.493 is Reff at n=1.4
x0 = z0.^2 + rho^2;
x1 = (z0+2.*zb).^2 + rho^2; 
model_refl = (1./4*pi).*((z0.*(mue+1./sqrt(x0)).*exp(-mue.*sqrt(x0))./x0) + ((z0 + 2*zb).*(mue + 1./sqrt(x1)).*exp(-mue.*sqrt(x1))./x1));


(pi*((exp(-(3*ext*hb*log10*(a/lambda^b + ext*hb*log10)*(1/(a/lambda^b + ext*hb*log10)^2 + rho^2)^(1/2))/2)*(1/(1/(a/lambda^b + ext*hb*log10)^2 + rho^2)^(1/2) + (3*ext*hb*log10*(a/lambda^b + ext*hb*log10))/2))/((a/lambda^b + ext*hb*log10)*(1/(a/lambda^b + ext*hb*log10)^2 + rho^2)) + (exp(-(3*ext*hb*log10*(a/lambda^b + ext*hb*log10)*((1/(a/lambda^b + ext*hb*log10) + 5972/(507*((3*a)/lambda^b + 3*ext*hb*log10)))^2 + rho^2)^(1/2))/2)*(1/(a/lambda^b + ext*hb*log10) + 5972/(507*((3*a)/lambda^b + 3*ext*hb*log10)))*(1/((1/(a/lambda^b + ext*hb*log10) + 5972/(507*((3*a)/lambda^b + 3*ext*hb*log10)))^2 + rho^2)^(1/2) + (3*ext*hb*log10*(a/lambda^b + ext*hb*log10))/2))/((1/(a/lambda^b + ext*hb*log10) + 5972/(507*((3*a)/lambda^b + 3*ext*hb*log10)))^2 + rho^2)))/4




for i = 1:length(Time)
    y = raw(i,:);

    Boundaries = [0.95,60,60,40,1;0.7,10,10,1,0.1];
    % Boundaries = [0.85, 24, 60, 40, 0.8; 0.70, 10, 10, 25, 0.3]; % this gives physiological upper and lower boundaries of WF, HHb, HbO2
    % Water Fraction in %, HHb and HbO2 in uM

    %WF fitting
    wave_start = 825; %starting wavelength
    wave_end = 850;
    exp_refl = generateexp_refl(wl_cal, x, y, wave_start, wave_end); %generate experimental reflectance in the 825-850 region
    diff_2_exp_refl = diff(smooth(diff(smooth(exp_refl)))); %take 2nd diff of exp refl
    xdata = wave_start:wave_end;

    lb = Boundaries(2,:);
    ub = Boundaries(1,:);
    % here xdata corresponds to the wavelength range we want to use
    diff_2_model_refl = @(parameters,xdata) gen_refl_model_diff2_WF(parameters ,xdata); % Second derivative of Model refl
    options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-25,'FunctionTolerance',1e-25,'StepTolerance',1e-20,'MaxFunctionEvaluations',100000,'MaxIterations',10000,'Display','off');
    [fit_parameters_1,resnorm] = lsqcurvefit(diff_2_model_refl,[0.60 10 20 30 0.6], xdata, diff_2_exp_refl, lb, ub,options);
    WF = fit_parameters_1(1);
    % X = LSQCURVEFIT(FUN(X,XDATA),X0,XDATA,YDATA)
    % parameters = [WF, HHb, HbO2, a, b]

    % HHb fitting
    wave_start = 700;
    wave_end = 800;
    exp_refl = generateexp_refl(wl_cal, x, y, wave_start, wave_end);
    diff_2_exp_refl = diff(smooth(diff(smooth(exp_refl)))); %take 2nd diff of exp refl
    xdata = wave_start:wave_end;

    lb = Boundaries(2,:);
    ub = Boundaries(1,:);
    lb(1,1) = WF;
    ub(1,1) = WF;
    diff_2_model_refl = @(parameters,xdata) gen_refl_model_diff2_HHb(parameters ,xdata, WF); % first derivative of model refl
    options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-25,'FunctionTolerance',1e-25,'StepTolerance',1e-20,'MaxFunctionEvaluations',10000,'MaxIterations',1000,'Display','off');
    fit_parameters_2 = lsqcurvefit(diff_2_model_refl,[0.60 10 20 30 0.6], xdata, diff_2_exp_refl, lb, ub, options);
    HHb = fit_parameters_2(2);
    xdata_2nd = xdata;

    %HbO2 fitting
    wave_start = 680;
    wave_end = 850;
    exp_refl = generateexp_refl(wl_cal, x, y, wave_start, wave_end);
    diff_1_exp_refl = diff(smooth(exp_refl)); %take 1st diff of exp refl
    lb(1,2) = HHb;
    ub(1,2) = HHb;
    xdata = wave_start:wave_end;
    diff_1_model_refl = @(parameters,xdata) gen_refl_model_diff1_HbO2(parameters,xdata,WF,HHb); % first derivative of model refl
    options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-25,'FunctionTolerance',1e-25,'StepTolerance',1e-20,'MaxFunctionEvaluations',10000,'MaxIterations',1000,'Display','off');
    fit_parameters_3 = lsqcurvefit(diff_1_model_refl,[0.60 10 20 30 0.6], xdata, diff_1_exp_refl, lb, ub, options);
    HbO2 = fit_parameters_3(3);
    a = fit_parameters_3(4); % physiological boundaries of a and b need to be added
    b = fit_parameters_3(5);


    Abs.WF(i) = WF;
    Abs.HHb(i) = HHb;
    Abs.HbO2(i) = HbO2;
    Abs.a(i) = a;
    Abs.b(i) = b;

    StO2(i) = HbO2/(HbO2+HHb)*100;

end

%% calculate dynamic pathlength

if PL_method == 2
    DPF_water = (Abs.WF./WaterFraction)*10./Settings.optode_dist;
    PL_water = DPF_water.*Settings.optode_dist;
else
    PL_water = [];
end

%% scale Conc to pathlength
if PL_method == 2
    Conc_noPL = Conc*(Settings.optode_dist*Settings.DPF);
    Conc_PL_water = Conc_noPL.*PL_water';
end
    