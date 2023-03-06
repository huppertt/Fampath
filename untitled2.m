 methods (Access = public)
    
    function [conc] = UCLn_code_NeoLight_GUI(app, raw, raw0, wl_cal, ext_coeffs_inv, pl, wl)
        
    %function to run UCLn algorithm to calculate concentrations (Bale et al. 2016 CCO Review)

    %set pathlength     
    if pl == 1
        lambda_dep = app.epsilon_labview_770to906(10:130,5);
    elseif pl ==0
        lambda_dep = ones(121,1);
    end


    %calculate attenuation
    atten = log10(raw0./raw); % log 10
    atten_int(:,1) = interp1(wl_cal, atten(1,:), wl, 'spline'); % interpolate to wavelengths of interest
    atten_int_wldep = atten_int ./ lambda_dep; % take into account wavelength dependency

    %calculate concentration from modified Beer Lambert Law
    conc =  double(ext_coeffs_inv) * double(atten_int_wldep) .* 1/(double(app.optode_dist)*double(app.DPF));
        
    end
        
    end
    
    
    
        methods (Access = public)
        
    function [WF, HHb, HbO2, a, b] = BF_function(app, wavelengths, ref, data)
    %function to calculate tissue saturation from broadband fitting of derivatives
    
    warning('off','all')
    
    app.UIAxes.cla;

    %load data from miniCYRIL 
    % wavelengths = Data.Wavelengths; 
    x = app.ref; 
    y = data; 
    
    Boundaries = [0.95,60,60,40,1;0.7,10,10,1,0.1];
    % Boundaries = [0.85, 24, 60, 40, 0.8; 0.70, 10, 10, 25, 0.3]; % this gives physiological upper and lower boundaries of WF, HHb, HbO2
    % Water Fraction in %, HHb and HbO2 in uM
    
    %WF fitting 
    wave_start = 825; %starting wavelength
    wave_end = 850;
    wl = double(wavelengths);
    exp_refl = generateexp_refl(wl, x, y, wave_start, wave_end); %generate experimental reflectance in the 825-850 region
    diff_2_exp_refl = diff(smooth(diff(smooth(exp_refl)))); %take 2nd diff of exp refl
    xdata = wave_start:wave_end;
    
    lb = Boundaries(2,:); 
    ub = Boundaries(1,:);
    % here xdata corresponds to the wavelength range we want to use 
    diff_2_model_refl = @(parameters,xdata) gen_refl_model_diff2_WF(parameters ,xdata); % Second derivative of Model refl 
    options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-25,'FunctionTolerance',1e-25,'StepTolerance',1e-20,'MaxFunctionEvaluations',100000,'MaxIterations',10000, 'Display', 'off');
    [fit_parameters_1,resnorm] = lsqcurvefit(diff_2_model_refl,[0.60 10 20 30 0.6], xdata, diff_2_exp_refl, lb, ub,options);
    WF = fit_parameters_1(1);
    % X = LSQCURVEFIT(FUN(X,XDATA),X0,XDATA,YDATA)
    % parameters = [WF, HHb, HbO2, a, b]
    
    % HHb fitting 
    wave_start = 700; 
    wave_end = 800; 
    exp_refl = generateexp_refl(wl, x, y, wave_start, wave_end); 
    diff_2_exp_refl = diff(smooth(diff(smooth(exp_refl)))); %take 2nd diff of exp refl
    xdata = wave_start:wave_end;
    
    lb = Boundaries(2,:); 
    ub = Boundaries(1,:);
    lb(1,1) = WF;
    ub(1,1) = WF;
    diff_2_model_refl = @(parameters,xdata) gen_refl_model_diff2_HHb(parameters ,xdata, WF); % first derivative of model refl 
    options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-25,'FunctionTolerance',1e-25,'StepTolerance',1e-20,'MaxFunctionEvaluations',10000,'MaxIterations',1000, 'Display', 'off');
    fit_parameters_2 = lsqcurvefit(diff_2_model_refl,[0.60 10 20 30 0.6], xdata, diff_2_exp_refl, lb, ub, options); 
    HHb = fit_parameters_2(2); 
    xdata_2nd = xdata; 
    
    %HbO2 fitting 
    wave_start = 680; 
    wave_end = 850; 
    exp_refl = generateexp_refl(wl, x, y, wave_start, wave_end); 
    diff_1_exp_refl = diff(smooth(exp_refl)); %take 1st diff of exp refl
    lb(1,2) = HHb; 
    ub(1,2) = HHb;
    xdata = wave_start:wave_end;
    diff_1_model_refl = @(parameters,xdata) gen_refl_model_diff1_HbO2(parameters,xdata,WF,HHb); % first derivative of model refl 
    options = optimoptions('lsqcurvefit','OptimalityTolerance',1e-25,'FunctionTolerance',1e-25,'StepTolerance',1e-20,'MaxFunctionEvaluations',10000,'MaxIterations',1000, 'Display', 'off');
    fit_parameters_3 = lsqcurvefit(diff_1_model_refl,[0.60 10 20 30 0.6], xdata, diff_1_exp_refl, lb, ub, options); 
    HbO2 = fit_parameters_3(3); 
    a = fit_parameters_3(4); % physiological boundaries of a and b need to be added 
    b = fit_parameters_3(5); 

    end

    end
    
    
    
  
            
            
