%% Calculate concentration from back up miniCYRIL data

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters to check before running script%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PL_method = 2;
% 1 = DPF, 2 = water fitting

%for DPF
optode_dist = 3;
DPF = 4.99; %adult forearm 4.16, baby head 4.99, adult head 6.26, adult leg 5.51 (Duncan et al. 1995)

%for water fitting
WaterFraction = 0.85;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data

[file, path] = uigetfile('C:\Users\Justin\Documents\miniCYRIL in the Antipodes\Back Up NIRS Data\*.txt','Select back up file');
s_and_time = load([path file]);

filename = ['C:\Users\Justin\Documents\miniCYRIL in the Antipodes\Processed NIRS Data\' erase(file,'.txt') ' processed.mat'];

Spectra = s_and_time(2:end,2:end);
Time = s_and_time(2:end,1);

raw = Spectra; 

%% Start calculation
% Read in raw CYRIL files
% to be written as necessary %

% Read in calibration file
wl_cal = load('wavelengths_Aus.mat');
wl_cal = wl_cal.Wavelengths;
Wavelengths = wl_cal;
wl = 780:900; % Wavelengths used to resolve over

% Extinction coefficients
c=tissue_specific_extinction_coefficient_650to1000;

% HbO2, HHb, CCO
ext_coeffs = c(wl(1)-649:wl(end)-649,3:5);
ext_coeffs_inv = pinv(ext_coeffs);

DPF_dep=DPF_Lambda_Dependency_740to915;
DPF_dep = DPF_dep(wl(1)-739:wl(end)-739,2);

    
    %% calculate tissue saturation
    warning('off','all')

    ref1 = load('ref.mat');
    x = ref1.ref_save;
     
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
    DPF_water = (Abs.WF./WaterFraction)*10./optode_dist;
    PL_water = DPF_water.*optode_dist;
    else
        PL_water = [];
    end
    
    %% Calculate concentration
    
    % attenuation

[x,y] = size(raw);

for ii = 1:x % row
    for jj = 1:y % column
        atten(ii,jj) = log10(raw(1,jj)/raw(ii,jj)); % log 10
    end

atten_int(:,ii) = interp1(wl_cal, atten(ii,:), wl, 'spline'); % interpolate to wavelengths of interest

end

for ii = 1:length(wl)
    atten_int_wldep(ii,:) = atten_int(ii,:) / DPF_dep(ii); % take into account wavelength dependency
end

% Calculate concentration
for ii = 1:x
    if PL_method == 1
    Conc(ii,:) =  ext_coeffs_inv * atten_int_wldep(:,ii) .* 1/(optode_dist*DPF);
    else
    Conc(ii,:) =  ext_coeffs_inv * atten_int_wldep(:,ii) .* 1/(PL_water(ii));  
    end
end
    
    for i = 2:length(Time)
        Diff(i) = Time(i)-Time(i-1);
    end
    
    %% save

    Settings.DPF = DPF;
    Settings.optode_dist = optode_dist;
    Setting.PL_method = PL_method;

    clearvars -except PL_water Conc Wavelengths Time Spectra Settings Diff filename Abs StO2
    
    save(filename,'Spectra','Conc','Settings','Wavelengths','Time','Diff','Abs','StO2','PL_water')
