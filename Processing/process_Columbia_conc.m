function [HbO2,Hb,CCO,Abs]=process_Columbia_conc(Spectra,wavelengths,Ref,PL_method)

%Spatial-temporal PCA filtering
[U,S,V]=nirs.math.mysvd(Spectra);
lst=nirs.math.sig_eigens(diag(S));
Spectra=U(:,lst)*S(lst,lst)*V(:,lst)';

if(nargin<4)
    PL_method = 1;
    % 1 = DPF, 2 = water fitting
end

Settings.optode_dist=30;
Settings.DPF=6.26;

%for water fitting
WaterFraction = 0.85;

wl = 650:915;
Ref = interp1(wavelengths,Ref,wl,'linear');

 Spectra2=[];
for i=1:size(Spectra,1)
    Spectra2(i,:) = interp1(wavelengths,Spectra(i,:),wl,'linear');
end
Spectra=Spectra2;



c=tissue_specific_extinction_coefficient_650to1000;
% HbO2, HHb, CCO
ext_coeffs=[];
for i=3:5
    ext_coeffs(:,i-2) = interp1(c(:,1),c(:,i),wl);
end


DPF_dep=DPF_Lambda_Dependency_740to915;
DPF_dep = interp1(DPF_dep(:,1),DPF_dep(:,2),wl);
DPF_dep(isnan(DPF_dep))=1.0273;


     
Boundaries = [0.95,60,60,40,1;0.7,10,10,1,0.1];
% Boundaries = [0.85, 24, 60, 40, 0.8; 0.70, 10, 10, 25, 0.3]; % this gives physiological upper and lower boundaries of WF, HHb, HbO2
% Water Fraction in %, HHb and HbO2 in uM

Spectra_smooth=resample(Spectra,1,30);
Abs=struct;

if(PL_method==2)
%% calculate tissue saturation
for i = 1:size(Spectra_smooth,1)
    %WF fitting
    wave_start = 825; %starting wavelength
    wave_end = 850;
    exp_refl = generateexp_refl(wl, Ref, Spectra_smooth(i,:), wave_start, wave_end); %generate experimental reflectance in the 825-850 region
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
    exp_refl = generateexp_refl(wl, Ref, Spectra(i,:), wave_start, wave_end);
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
    exp_refl = generateexp_refl(wl, Ref, Spectra(i,:), wave_start, wave_end);
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


    Abs.WF(i,1) = WF;
    Abs.HHb(i,1) = HHb;
    Abs.HbO2(i,1) = HbO2;
    Abs.a(i,1) = a;
    Abs.b(i,1) = b;
    Abs.StO2(i,1) = HbO2/(HbO2+HHb)*100;
end
end
%% calculate dynamic pathlength

if PL_method == 2
    DPF_water = (Abs.WF./WaterFraction)./Settings.optode_dist;
    Abs.PL_water = (DPF_water.*Settings.optode_dist);
else
    PL_water = [];
end


%% Calculate concentration
atten = log10(Spectra./(ones(size(Spectra,1),1)*median(Spectra,1)));
atten_int_wldep = atten./(ones(size(Spectra,1),1)*DPF_dep);


ext_coeffs(:,4)=Ref;
ext_coeffs(:,5)=1;

w=1./var(atten_int_wldep);

ext_coeffs_inv = pinv(ext_coeffs);
ext_coeffs_inv = pinv(ext_coeffs'*diag(w)*ext_coeffs)*ext_coeffs'*diag(w);


% Calculate concentration
if PL_method == 1
    Conc =  (ext_coeffs_inv * atten_int_wldep' /(Settings.optode_dist*Settings.DPF))'*1E6;
else
    PL_water=interp1(linspace(0,1,length(Abs.PL_water)),Abs.PL_water,linspace(0,1,size(atten,1)));
    Conc=  ((ext_coeffs_inv * atten_int_wldep')./(ones(size(ext_coeffs,2),1)*PL_water))'*1E5;
end

HbO2=Conc(:,1);
Hb=Conc(:,2);
CCO=Conc(:,3);



