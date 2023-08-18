classdef miniCYRIL_at_Columbia_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ReadyLampLabel                  matlab.ui.control.Label
        ReadyLamp                       matlab.ui.control.Lamp
        TabGroup                        matlab.ui.container.TabGroup
        SettingsTab                     matlab.ui.container.Tab
        GeneralPanel                    matlab.ui.container.Panel
        IntegrationTimesEditFieldLabel  matlab.ui.control.Label
        IntegrationTimesEditField       matlab.ui.control.NumericEditField
        GroupnameEditFieldLabel         matlab.ui.control.Label
        GroupnameEditField              matlab.ui.control.EditField
        PatientIDEditFieldLabel         matlab.ui.control.Label
        PatientIDEditField              matlab.ui.control.EditField
        SynchronisationDropDownLabel    matlab.ui.control.Label
        SynchronisationDropDown         matlab.ui.control.DropDown
        PathlengthPanel                 matlab.ui.container.Panel
        PathlengthMethodButtonGroup     matlab.ui.container.ButtonGroup
        DPFButton                       matlab.ui.control.RadioButton
        NopathlengthButton              matlab.ui.control.RadioButton
        WaterfittingButton              matlab.ui.control.RadioButton
        SeparationcmEditFieldLabel      matlab.ui.control.Label
        SeparationcmEditField           matlab.ui.control.NumericEditField
        DPFDropDownLabel                matlab.ui.control.Label
        DPFDropDown                     matlab.ui.control.DropDown
        WaterFractionEditFieldLabel     matlab.ui.control.Label
        WaterFractionEditField          matlab.ui.control.NumericEditField
        TemperaturePanel                matlab.ui.container.Panel
        SettempButton                   matlab.ui.control.Button
        SetTemperatureCEditFieldLabel   matlab.ui.control.Label
        SetTemperatureCEditField        matlab.ui.control.NumericEditField
        ChecktempButton                 matlab.ui.control.Button
        PathlengthcmEditFieldLabel      matlab.ui.control.Label
        PathlengthcmEditField           matlab.ui.control.NumericEditField
        ConcentrationsTab               matlab.ui.container.Tab
        UIAxes2                         matlab.ui.control.UIAxes
        NumbertodisplayEditFieldLabel   matlab.ui.control.Label
        NumbertodisplayEditField        matlab.ui.control.NumericEditField
        SpectraTab                      matlab.ui.container.Tab
        UIAxes                          matlab.ui.control.UIAxes
        ReferenceTab                    matlab.ui.container.Tab
        UIAxes_2                        matlab.ui.control.UIAxes
        ReferencespectrawascollectedonTextAreaLabel  matlab.ui.control.Label
        ReferencespectrawascollectedonTextArea  matlab.ui.control.TextArea
        ImportrefButton                 matlab.ui.control.Button
        TissueSaturationTab             matlab.ui.container.Tab
        UIAxes3                         matlab.ui.control.UIAxes
        InitialiseButton                matlab.ui.control.Button
        nEditFieldLabel                 matlab.ui.control.Label
        nEditField                      matlab.ui.control.NumericEditField
        ControlsPanel                   matlab.ui.container.Panel
        StartButton                     matlab.ui.control.Button
        MarkEventButton                 matlab.ui.control.Button
        StopButton                      matlab.ui.control.Button
        EventDetailsEditFieldLabel      matlab.ui.control.Label
        EventDetailsEditField           matlab.ui.control.EditField
        LogeventdetailsButton           matlab.ui.control.Button
        miniCYRILatColumbiaLabel        matlab.ui.control.Label
        LiveTemperatureCEditFieldLabel  matlab.ui.control.Label
        LiveTemperatureCEditField       matlab.ui.control.NumericEditField
        LiveTimesEditFieldLabel         matlab.ui.control.Label
        LiveTimesEditField              matlab.ui.control.NumericEditField
    end


    properties (Access = private)
       
        % load the drivers for Wasatch spectrometer
        dll = NET.addAssembly('C:\Program Files\Wasatch Photonics\Wasatch.NET\WasatchNET.dll');
        driver = WasatchNET.Driver.getInstance();
        
        %set defaults
        stopcommand = 0;
        event_button = 0;
        initialised = 0;
        n = 36000; %preallocation number
        Spectrums = [];
        SpectrumInit = [];
        DPF = [];
        optode_dist = [];
        StO2 = [];
        ref = [];       
        r = 0;
        Abs = [];
        s = [];
        t = timer;
        e = [];
        event_str = [];
        % Extinction coefficients: Wavelength, HbO2, HHb, CCO, wavelength dependency
        epsilon_labview_770to906 = specific_extinction_coeffs_770to906;
 
    end 

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


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function getStarted(app)
        %function that sets defaults once app is opened
        
            %display reference
            ref1 = load('ref.mat');
            app.ref = ref1.ref_save;
            date = ref1.ts(1,:);
            wave1 = load('wavelengths_Col.mat');
            wavelengths = wave1.Wavelengths;
            r1 = line(app.UIAxes_2, wavelengths, app.ref);
            xlim(app.UIAxes_2, [600 950])
            ylim(app.UIAxes_2, [0 70000])
            
            app.ReferencespectrawascollectedonTextArea.Value = date;
            
            if datenum(date) <= (now - 30)
                uialert(app.UIFigure,'Reference spectra is more than a month old. Take new reference spectra and import ASAP.','Reference Spectra Warning')
             end
        
            %open spectrometer and set temperature
            numberOfSpectrometers = app.driver.openAllSpectrometers();

            if numberOfSpectrometers <= 0
            	uialert(app.UIFigure,'No spectrometer found. Check USB connection.','Spectrometer Error')
                app.ReadyLamp.Color = [1, 0, 0];
                app.ReadyLampLabel.Text = 'Not connected';
             return
            end
            
            spectrometer = app.driver.getSpectrometer(0);
            
            spectrometer.detectorTECEnabled = 1;
            spectrometer.detectorTECSetpointDegC = app.SetTemperatureCEditField.Value;                      
                       
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
        %main function when start button pushed - collects data and plots it

   app.ReadyLamp.Color = [0, 0, 1]; %set ready button to blue (running)
   app.ReadyLampLabel.Text = 'Collecting data';

   app.UIAxes2.cla;
   app.UIAxes.cla;
   app.UIAxes3.cla;

   app.stopcommand = 0; 
   
   if app.initialised == 0 %Warning: if not initialised and...
      if app.WaterfittingButton.Value ==1 %water fitting pathlengths selected 
             uialert(app.UIFigure,'Cannot calculate pathlength from water fitting because initialisation has not been performed. Initialise or choose different pathlength option.','Software Error')
             app.ReadyLamp.Color = [0, 1, 0];
             app.ReadyLampLabel.Text = 'Ready';
             return
      end
   end
 
   cla(app.UIAxes3)

    %open spectrometer
    numberOfSpectrometers = app.driver.openAllSpectrometers();
    
    if numberOfSpectrometers <= 0
       uialert(app.UIFigure,'No spectrometer found. Check USB connection.','Spectrometer Error')
       app.ReadyLamp.Color = [1, 0, 0];
       app.ReadyLampLabel.Text = 'Not connected';
       return
    end
    
    spectrometer = app.driver.getSpectrometer(0);
    spectrometer.detectorTECEnabled = 1;
    spectrometer.detectorTECSetpointDegC = app.SetTemperatureCEditField.Value; 

    %check temperature
    app.LiveTemperatureCEditField.Value = double(spectrometer.detectorTemperatureDegC); %display temp

    if app.LiveTemperatureCEditField.Value >= -14
    
        temp_selection = uiconfirm(app.UIFigure,'Be aware that the detector temperature is too high. It is recommended that you stop and wait until temperature reaches -15C. Do you want to continue?','Detector Temperature Warning',...
                        'Icon','warning','Options',{'Continue','Stop'},...
                        'DefaultOption',2,'CancelOption',2);
            
           if contains(temp_selection,'Stop')==1
                app.ReadyLamp.Color = [0, 1, 0];
                app.ReadyLampLabel.Text = 'Ready';
                return        
           end
     end

% access some key properties
pixels = spectrometer.pixels;
modelName = char(spectrometer.model);
serialNum = char(spectrometer.serialNumber);
wavelengths = spectrometer.wavelengths;

%set temperature
spectrometer.detectorTECEnabled = 1;
spectrometer.detectorTECSetpointDegC = app.SetTemperatureCEditField.Value;
spectrometer.integrationTimeMS = app.IntegrationTimesEditField.Value*1000;
Temperature = spectrometer.detectorTemperatureDegC;


%set up pathlength
%load DPF and optode distance (DPF from Duncan et al. 1995)
switch app.DPFDropDown.Value
                case 'Baby Head'
                app.DPF = 4.99;
                case 'Adult Head'
                app.DPF = 6.26;
                case 'Adult Arm'
                app.DPF = 4.16;
                case 'Adult Leg'   
                app.DPF = 5.51;
 end
        
app.optode_dist = app.SeparationcmEditField.Value;

if app.DPFButton.Value ==1
      
        DPF_code = app.DPF;
        optode_dist_code = app.optode_dist;
        pl =1;
        app.PathlengthcmEditField.Value = DPF_code*app.optode_dist;
        
elseif app.NopathlengthButton.Value ==1
    
    DPF_code = 1;
    optode_dist_code = 1;
    pl = 0;
    app.PathlengthcmEditField.Value = DPF_code*optode_dist_code;
    
elseif app.WaterfittingButton.Value ==1

    DPF_code = (app.Abs.WF/app.WaterFractionEditField.Value)*10/app.optode_dist;
    pl = 1;
    optode_dist_code =  app.optode_dist;
    app.PathlengthcmEditField.Value = DPF_code*app.optode_dist;
    
end

% get ext coeffs: HbO2, HHb, CCO
ext_coeffs = app.epsilon_labview_770to906(10:130,2:4);
ext_coeffs_inv = pinv(ext_coeffs);

wl = 780:900; %wavelengths to resolve over

overall = tic;
i = 1;
app.e = 1;

%set up back up file
DateString = datestr(datetime('now'));
DateString(regexp(DateString,':'))=[];
filename = ['C:\Users\uceegm0\Documents\miniCYRIL at Columbia\Raw NIRS Data\' DateString ' ' app.GroupnameEditField.Value ' Patient ' app.PatientIDEditField.Value];
filename_backup = ['C:\Users\uceegm0\Documents\miniCYRIL at Columbia\Back Up NIRS Data\' DateString ' ' app.GroupnameEditField.Value ' Patient ' app.PatientIDEditField.Value ' back up.txt'];

s_and_time=zeros(1,1025);
s_and_time(1,2:end) = wavelengths;

save(filename_backup,'s_and_time','-ascii');

%preallocate arrays
conc = nan(app.n,3);
ti = nan(app.n,1);
events = nan(app.n,1);
Temp = nan(app.n,1);
StO2_cont = nan(app.n,1);
app.Spectrums = nan(app.n,1024);
app.StO2 = nan(app.n,1);

% set up graphs

spectrum = nan(1024,1);

s1 = line(app.UIAxes, wavelengths, spectrum);
xlim(app.UIAxes, [600 950])
ylim(app.UIAxes, [0 70000])

 h1 = line(app.UIAxes2, 'XData', ti, 'YData', conc(:,1),'Color','r');
 h2 = line(app.UIAxes2, 'XData', ti, 'YData', conc(:,2),'Color','b');
 h3 = line(app.UIAxes2, 'XData', ti, 'YData', conc(:,3),'Color','g');
 legend(app.UIAxes2,[h1,h2,h3],{'HbO_2','HHb','oxCCO'},'Location','northwest')
 
 if app.initialised == 1
        st1 = line(app.UIAxes3,ti,StO2_cont);
 end

   if strcmp(app.SynchronisationDropDown.Value,'None')==0
       id = 'daq:Session:onDemandOnlyChannelsAdded';
        warning('off',id);
    
        d = daq.getDevices;
        s = daq.createSession('ni'); %or daq.createSession(d.Vendor.ID);
        ch = addDigitalChannel(s,'dev1','Port1/Line3','OutputOnly');
   end
 

while  app.stopcommand == 0 %main loop running to collect data and plot
 
       
    if i == 1
    if strcmp(app.SynchronisationDropDown.Value,'Start')==1
    
        
        outputSingleScan(s,1)
        release(s)

        
    end
    end
    
    if strcmp(app.SynchronisationDropDown.Value,'Every Measurement')==1
        
        outputSingleScan(s,1)
        release(s)
    end
    
    %get intensity spectrum
    spectrum = spectrometer.getSpectrum();
    
    if i == 1
    if strcmp(app.SynchronisationDropDown.Value,'Start')==1
    
        outputSingleScan(s,0)
        release(s)

        
    end
    end
    
    if strcmp(app.SynchronisationDropDown.Value,'Every Measurement')==1
        outputSingleScan(s,0)
        release(s)
    end

    ti(i) = toc(overall);
    Timestamp(i,:) = datestr(datetime('now'));
    
    %save back up file
    s_and_time=horzcat(toc(overall),double(spectrum));
    save(filename_backup,'s_and_time','-append','-ascii'); %make this append
 
    if ~isempty(spectrum)
        app.Spectrums(i,:)=double(spectrum);
        app.LiveTimesEditField.Value = ti(i);
    
        %check for events
        
        if app.event_button == 1
            events(i,1) = 1;
            events(i,2) = app.e;
            app.event_button = 0;
            app.e=app.e+1;
        else
            events(i,1:2) = 0;
        end
        
        %display temperature
        Temp(i,1) = double(spectrometer.detectorTemperatureDegC);
        app.LiveTemperatureCEditField.Value = double(spectrometer.detectorTemperatureDegC); %display temp
        
        %Calc conc
       
        if i == 1
        spectrum0 = double(spectrum);
        end
     
        conc(i,:) = UCLn_code_NeoLight_GUI(app, double(spectrum), spectrum0, double(wavelengths), ext_coeffs_inv, pl, wl);
      
        %plot data during cycle
        
        display = app.NumbertodisplayEditField.Value;
    
        if i<display+1
        set(h1, 'YData', conc(:,1));
        set(h1, 'XData', ti);
        set(h2, 'YData', conc(:,2));
        set(h2, 'XData', ti);
        set(h3, 'YData', conc(:,3));
        set(h3, 'XData', ti);
        else
        set(h1, 'YData', conc(i-display:i,1));
        set(h1, 'XData', ti(i-display:i));
        set(h2, 'YData', conc(i-display:i,2));
        set(h2, 'XData', ti(i-display:i));
        set(h3, 'YData', conc(i-display:i,3));
        set(h3, 'XData', ti(i-display:i));
        end
        
        set(s1, 'YData', spectrum);
              
        if app.initialised == 1
     
            StO2_cont(i,1) = (app.Abs.HbO2 + conc(i,1)*1000)/(app.Abs.HbO2 + conc(i,1)*1000 + app.Abs.HHb + conc(i,2)*1000)*100;
    
            if i<display+1
            set(st1,'XData',ti);
            set(st1,'YData',StO2_cont);
            else
            set(st1,'XData',ti(i-display:i));
            set(st1,'YData',StO2_cont(i-display:i));
            end
            
        end
        
        pause(0.1)
        
        i = i+1;
     
    end
end

%save Settings
Settings.DPF = DPF_code;
Settings.optode_dist = optode_dist_code;
Settings.integration_time_ms = app.IntegrationTimesEditField.Value*1000;

%save data

Conc = rmmissing(conc);
Abs = app.Abs;
StO2 = rmmissing(StO2_cont);
Wavelengths = double(wavelengths);
Time = rmmissing(ti);
Events = rmmissing(events);
Event_details = app.event_str;
Spectra = rmmissing(app.Spectrums);
Ref = app.ref;
Temp = rmmissing(Temp);

for i = 2:length(Time)
    Diff(i) = Time(i)-Time(i-1);
end

save(filename,'Spectra','Conc','Temp','Settings','Abs','StO2','Event_details','Wavelengths','Time','Events','Diff','Ref','Timestamp')
                   
clear Conc Abs StO2 Event_details Wavelengths Time Spectra Temp Settings Events Diff Ref

app.ReadyLamp.Color = [0, 1, 0];
app.ReadyLampLabel.Text = 'Ready';

        end

        % Value changed function: IntegrationTimesEditField
        function IntegrationTimesEditFieldValueChanged(app, event)
            spectrometer.integrationTimeMS = app.IntegrationTimesEditField.Value*1000;
        end

        % Value changed function: SetTemperatureCEditField
        function SetTemperatureCEditFieldValueChanged(app, event)
            spectrometer.detectorTECSetpointDegC = app.SetTemperatureCEditField.Value;
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
                   app.stopcommand = 1; 
        end

        % Button pushed function: InitialiseButton
        function InitialiseButtonPushed(app, event)
        %acquire n seconds of spectrum and apply broadband fitting to get StO2 
        
            app.ReadyLamp.Color = [0, 0, 1];
            app.ReadyLampLabel.Text = 'Collecting data';
            
            %open spectrometer
            
            numberOfSpectrometers = app.driver.openAllSpectrometers();

            if numberOfSpectrometers <= 0
                   uialert(app.UIFigure,'No spectrometer found. Check USB connection.','Spectrometer Error')
                   return
            end
            
            spectrometer = app.driver.getSpectrometer(0);

            %check temperature
            app.LiveTemperatureCEditField.Value = double(spectrometer.detectorTemperatureDegC); %display temp

   
            if app.LiveTemperatureCEditField.Value >= -14       
                temp_selection = uiconfirm(app.UIFigure,'Be aware that the detector temperature is too high. It is recommended that you stop and wait until temperature reaches -15C. Do you want to continue?','Detector Temperature Warning',...
                                        'Icon','warning','Options',{'Continue','Stop'},...
                           'DefaultOption',2,'CancelOption',2);
                if contains(temp_selection,'Stop')==1
                    
                    app.ReadyLamp.Color = [0, 1, 0];
                    app.ReadyLampLabel.Text = 'Ready';
                    return
                end
            end
 
% access some key properties
pixels = spectrometer.pixels;
modelName = char(spectrometer.model);
serialNum = char(spectrometer.serialNumber);
wavelengths = spectrometer.wavelengths;

%set temperature
spectrometer.detectorTECEnabled = 1;
spectrometer.detectorTECSetpointDegC = app.SetTemperatureCEditField.Value;
spectrometer.integrationTimeMS = app.IntegrationTimesEditField.Value*1000;
Temperature = spectrometer.detectorTemperatureDegC;

DateString = datestr(datetime('now'));

tic      
for j = 1:app.nEditField.Value %loop to collect initial spectra 
    
    pause(0.00001)

    spectrumInit = spectrometer.getSpectrum();
    s_and_time=horzcat(toc,double(spectrumInit));
    app.SpectrumInit(j,:)=s_and_time;
    
    TempInit(j,1) = double(spectrometer.detectorTemperatureDegC);
    app.LiveTemperatureCEditField.Value = TempInit(j); %display temp
    pause(0.0001)
    
    %plot data during cycle
    line(app.UIAxes, wavelengths, spectrumInit);
    
    pause(0.00001)

end

    %broadband fitting    
    data = mean(app.SpectrumInit(:,2:end),1);
    
    if app.r == 0
        ref1 = load('ref.mat');
        app.ref = ref1.ref_save;
    end
                  
    [WF, HHb, HbO2, a, b] = BF_function(app, wavelengths, app.ref, data);
    
    %calculate tissue saturation
    app.StO2 = HbO2/(HbO2+HHb)*100;
    
    %plot initial StO2
    plot(app.UIAxes3,0,app.StO2,'*');
    
    app.Abs.WF = WF;
    app.Abs.HHb = HHb;
    app.Abs.HbO2 = HbO2;
    app.Abs.a = a;
    app.Abs.b = b;
    
    app.initialised = 1;
    app.ReadyLamp.Color = [0, 1, 0];
    app.ReadyLampLabel.Text = 'Ready';
 
        end

        % Button pushed function: ChecktempButton
        function ChecktempButtonPushed(app, event)
numberOfSpectrometers = app.driver.openAllSpectrometers();

if numberOfSpectrometers <= 0
	return
end

spectrometer = app.driver.getSpectrometer(0);

app.LiveTemperatureCEditField.Value = double(spectrometer.detectorTemperatureDegC); %display temp
            
        end

        % Button pushed function: SettempButton
        function SettempButtonPushed(app, event)
            
numberOfSpectrometers = app.driver.openAllSpectrometers();

if numberOfSpectrometers <= 0
	return
end

% open the first spectrometer found
spectrometer = app.driver.getSpectrometer(0);

%set temperature
spectrometer.detectorTECEnabled = 1;
spectrometer.detectorTECSetpointDegC = app.SetTemperatureCEditField.Value;
            
        end

        % Button pushed function: ImportrefButton
        function ImportrefButtonPushed(app, event)
ref_file = uigetfile();
app.UIFigure.Visible = 'off';
app.UIFigure.Visible = 'on';
ref_data = load(ref_file);
app.ref = nanmean(ref_data.Spectra,1);
ts = ref_data.Timestamp;
ref_save = app.ref;
save('ref.mat','ref_save','ts');
app.r = 1;
        end

        % Button pushed function: MarkEventButton
        function MarkEventButtonPushed(app, event)
app.event_button = 1;
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            
            %warning to check if want to quit
            close_selection = uiconfirm(app.UIFigure,'Are you sure you want to quit?','Close Warning',...
                        'Icon','warning','Options',{'Quit','Cancel'},...
           'DefaultOption',2,'CancelOption',2);
            
           if contains(close_selection,'Quit')
            delete(app)
           end         
        end

        % Button pushed function: LogeventdetailsButton
        function LogeventdetailsButtonPushed(app, event)
            % log details for previous event        
            if app.e>1
            app.event_str{app.e-1} = app.EventDetailsEditField.Value;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 958 480];
            app.UIFigure.Name = 'UI Figure';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create ReadyLampLabel
            app.ReadyLampLabel = uilabel(app.UIFigure);
            app.ReadyLampLabel.HorizontalAlignment = 'right';
            app.ReadyLampLabel.Position = [748 439 144 22];
            app.ReadyLampLabel.Text = 'Ready';

            % Create ReadyLamp
            app.ReadyLamp = uilamp(app.UIFigure);
            app.ReadyLamp.Position = [907 440 21 21];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [150 53 788 374];

            % Create SettingsTab
            app.SettingsTab = uitab(app.TabGroup);
            app.SettingsTab.Title = 'Settings';

            % Create GeneralPanel
            app.GeneralPanel = uipanel(app.SettingsTab);
            app.GeneralPanel.Title = 'General';
            app.GeneralPanel.Position = [23 32 260 276];

            % Create IntegrationTimesEditFieldLabel
            app.IntegrationTimesEditFieldLabel = uilabel(app.GeneralPanel);
            app.IntegrationTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.IntegrationTimesEditFieldLabel.Position = [26 216 109 22];
            app.IntegrationTimesEditFieldLabel.Text = 'Integration Time (s)';

            % Create IntegrationTimesEditField
            app.IntegrationTimesEditField = uieditfield(app.GeneralPanel, 'numeric');
            app.IntegrationTimesEditField.LowerLimitInclusive = 'off';
            app.IntegrationTimesEditField.Limits = [0 120];
            app.IntegrationTimesEditField.ValueChangedFcn = createCallbackFcn(app, @IntegrationTimesEditFieldValueChanged, true);
            app.IntegrationTimesEditField.Position = [150 216 100 22];
            app.IntegrationTimesEditField.Value = 1;

            % Create GroupnameEditFieldLabel
            app.GroupnameEditFieldLabel = uilabel(app.GeneralPanel);
            app.GroupnameEditFieldLabel.HorizontalAlignment = 'right';
            app.GroupnameEditFieldLabel.Position = [63 143 72 22];
            app.GroupnameEditFieldLabel.Text = 'Group name';

            % Create GroupnameEditField
            app.GroupnameEditField = uieditfield(app.GeneralPanel, 'text');
            app.GroupnameEditField.Position = [150 143 100 22];
            app.GroupnameEditField.Value = 'MIND';

            % Create PatientIDEditFieldLabel
            app.PatientIDEditFieldLabel = uilabel(app.GeneralPanel);
            app.PatientIDEditFieldLabel.HorizontalAlignment = 'right';
            app.PatientIDEditFieldLabel.Position = [77 115 58 22];
            app.PatientIDEditFieldLabel.Text = 'Patient ID';

            % Create PatientIDEditField
            app.PatientIDEditField = uieditfield(app.GeneralPanel, 'text');
            app.PatientIDEditField.Position = [150 115 100 22];
            app.PatientIDEditField.Value = '1';

            % Create SynchronisationDropDownLabel
            app.SynchronisationDropDownLabel = uilabel(app.GeneralPanel);
            app.SynchronisationDropDownLabel.HorizontalAlignment = 'right';
            app.SynchronisationDropDownLabel.Position = [33 48 91 22];
            app.SynchronisationDropDownLabel.Text = 'Synchronisation';

            % Create SynchronisationDropDown
            app.SynchronisationDropDown = uidropdown(app.GeneralPanel);
            app.SynchronisationDropDown.Items = {'Start', 'Every Measurement', 'None'};
            app.SynchronisationDropDown.Position = [139 48 100 22];
            app.SynchronisationDropDown.Value = 'Start';

            % Create PathlengthPanel
            app.PathlengthPanel = uipanel(app.SettingsTab);
            app.PathlengthPanel.Title = 'Pathlength';
            app.PathlengthPanel.Position = [295 65 198 243];

            % Create PathlengthMethodButtonGroup
            app.PathlengthMethodButtonGroup = uibuttongroup(app.PathlengthPanel);
            app.PathlengthMethodButtonGroup.Title = 'Pathlength Method';
            app.PathlengthMethodButtonGroup.Position = [67 124 123 91];

            % Create DPFButton
            app.DPFButton = uiradiobutton(app.PathlengthMethodButtonGroup);
            app.DPFButton.Text = 'DPF';
            app.DPFButton.Position = [11 45 58 22];
            app.DPFButton.Value = true;

            % Create NopathlengthButton
            app.NopathlengthButton = uiradiobutton(app.PathlengthMethodButtonGroup);
            app.NopathlengthButton.Text = 'No pathlength';
            app.NopathlengthButton.Position = [11 23 97 22];

            % Create WaterfittingButton
            app.WaterfittingButton = uiradiobutton(app.PathlengthMethodButtonGroup);
            app.WaterfittingButton.Text = 'Water fitting';
            app.WaterfittingButton.Position = [11 1 86 22];

            % Create SeparationcmEditFieldLabel
            app.SeparationcmEditFieldLabel = uilabel(app.PathlengthPanel);
            app.SeparationcmEditFieldLabel.HorizontalAlignment = 'right';
            app.SeparationcmEditFieldLabel.Position = [35 91 94 22];
            app.SeparationcmEditFieldLabel.Text = ' Separation (cm)';

            % Create SeparationcmEditField
            app.SeparationcmEditField = uieditfield(app.PathlengthPanel, 'numeric');
            app.SeparationcmEditField.Limits = [0.1 10];
            app.SeparationcmEditField.Position = [135 91 55 22];
            app.SeparationcmEditField.Value = 3;

            % Create DPFDropDownLabel
            app.DPFDropDownLabel = uilabel(app.PathlengthPanel);
            app.DPFDropDownLabel.HorizontalAlignment = 'right';
            app.DPFDropDownLabel.Position = [46 48 29 22];
            app.DPFDropDownLabel.Text = 'DPF';

            % Create DPFDropDown
            app.DPFDropDown = uidropdown(app.PathlengthPanel);
            app.DPFDropDown.Items = {'Baby Head', 'Adult Head', 'Adult Arm', 'Adult Leg'};
            app.DPFDropDown.Position = [90 48 100 22];
            app.DPFDropDown.Value = 'Baby Head';

            % Create WaterFractionEditFieldLabel
            app.WaterFractionEditFieldLabel = uilabel(app.PathlengthPanel);
            app.WaterFractionEditFieldLabel.HorizontalAlignment = 'right';
            app.WaterFractionEditFieldLabel.Position = [20 15 83 22];
            app.WaterFractionEditFieldLabel.Text = 'Water Fraction';

            % Create WaterFractionEditField
            app.WaterFractionEditField = uieditfield(app.PathlengthPanel, 'numeric');
            app.WaterFractionEditField.Limits = [0 1];
            app.WaterFractionEditField.Position = [118 15 70 22];
            app.WaterFractionEditField.Value = 0.85;

            % Create TemperaturePanel
            app.TemperaturePanel = uipanel(app.SettingsTab);
            app.TemperaturePanel.Title = 'Temperature';
            app.TemperaturePanel.Position = [510 32 260 276];

            % Create SettempButton
            app.SettempButton = uibutton(app.TemperaturePanel, 'push');
            app.SettempButton.ButtonPushedFcn = createCallbackFcn(app, @SettempButtonPushed, true);
            app.SettempButton.BackgroundColor = [0.8 0.8 0.8];
            app.SettempButton.Position = [136 169 100 22];
            app.SettempButton.Text = 'Set temp';

            % Create SetTemperatureCEditFieldLabel
            app.SetTemperatureCEditFieldLabel = uilabel(app.TemperaturePanel);
            app.SetTemperatureCEditFieldLabel.HorizontalAlignment = 'right';
            app.SetTemperatureCEditFieldLabel.Position = [7 216 114 22];
            app.SetTemperatureCEditFieldLabel.Text = 'Set Temperature (C)';

            % Create SetTemperatureCEditField
            app.SetTemperatureCEditField = uieditfield(app.TemperaturePanel, 'numeric');
            app.SetTemperatureCEditField.Limits = [-20 20];
            app.SetTemperatureCEditField.ValueChangedFcn = createCallbackFcn(app, @SetTemperatureCEditFieldValueChanged, true);
            app.SetTemperatureCEditField.Position = [136 216 100 22];
            app.SetTemperatureCEditField.Value = -15;

            % Create ChecktempButton
            app.ChecktempButton = uibutton(app.TemperaturePanel, 'push');
            app.ChecktempButton.ButtonPushedFcn = createCallbackFcn(app, @ChecktempButtonPushed, true);
            app.ChecktempButton.BackgroundColor = [0.902 0.902 0.902];
            app.ChecktempButton.Position = [14 169 100 22];
            app.ChecktempButton.Text = 'Check temp';

            % Create PathlengthcmEditFieldLabel
            app.PathlengthcmEditFieldLabel = uilabel(app.SettingsTab);
            app.PathlengthcmEditFieldLabel.HorizontalAlignment = 'right';
            app.PathlengthcmEditFieldLabel.Position = [307 32 90 22];
            app.PathlengthcmEditFieldLabel.Text = 'Pathlength (cm)';

            % Create PathlengthcmEditField
            app.PathlengthcmEditField = uieditfield(app.SettingsTab, 'numeric');
            app.PathlengthcmEditField.Position = [412 32 69 22];

            % Create ConcentrationsTab
            app.ConcentrationsTab = uitab(app.TabGroup);
            app.ConcentrationsTab.Title = 'Concentrations';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.ConcentrationsTab);
            title(app.UIAxes2, '')
            xlabel(app.UIAxes2, 'Time (s) ')
            ylabel(app.UIAxes2, 'Change (mM)')
            app.UIAxes2.PlotBoxAspectRatio = [1.96850393700787 1 1];
            app.UIAxes2.Position = [1 1 729 348];

            % Create NumbertodisplayEditFieldLabel
            app.NumbertodisplayEditFieldLabel = uilabel(app.ConcentrationsTab);
            app.NumbertodisplayEditFieldLabel.HorizontalAlignment = 'right';
            app.NumbertodisplayEditFieldLabel.Position = [662 205 102 22];
            app.NumbertodisplayEditFieldLabel.Text = 'Number to display';

            % Create NumbertodisplayEditField
            app.NumbertodisplayEditField = uieditfield(app.ConcentrationsTab, 'numeric');
            app.NumbertodisplayEditField.Limits = [1 Inf];
            app.NumbertodisplayEditField.Position = [661 182 100 22];
            app.NumbertodisplayEditField.Value = Inf;

            % Create SpectraTab
            app.SpectraTab = uitab(app.TabGroup);
            app.SpectraTab.Title = 'Spectra';

            % Create UIAxes
            app.UIAxes = uiaxes(app.SpectraTab);
            xlabel(app.UIAxes, 'Wavelength (nm)')
            ylabel(app.UIAxes, 'Intensity ')
            app.UIAxes.PlotBoxAspectRatio = [1.953125 1 1];
            app.UIAxes.Position = [1 1 788 348];

            % Create ReferenceTab
            app.ReferenceTab = uitab(app.TabGroup);
            app.ReferenceTab.Title = 'Reference';

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.ReferenceTab);
            xlabel(app.UIAxes_2, 'Wavelength (nm)')
            ylabel(app.UIAxes_2, 'Intensity ')
            app.UIAxes_2.PlotBoxAspectRatio = [1.953125 1 1];
            app.UIAxes_2.Position = [1 50 788 299];

            % Create ReferencespectrawascollectedonTextAreaLabel
            app.ReferencespectrawascollectedonTextAreaLabel = uilabel(app.ReferenceTab);
            app.ReferencespectrawascollectedonTextAreaLabel.HorizontalAlignment = 'right';
            app.ReferencespectrawascollectedonTextAreaLabel.Position = [155 11 200 22];
            app.ReferencespectrawascollectedonTextAreaLabel.Text = 'Reference spectra was collected on:';

            % Create ReferencespectrawascollectedonTextArea
            app.ReferencespectrawascollectedonTextArea = uitextarea(app.ReferenceTab);
            app.ReferencespectrawascollectedonTextArea.Position = [370 11 264 22];

            % Create ImportrefButton
            app.ImportrefButton = uibutton(app.ReferenceTab, 'push');
            app.ImportrefButton.ButtonPushedFcn = createCallbackFcn(app, @ImportrefButtonPushed, true);
            app.ImportrefButton.BackgroundColor = [0.8 0.8 0.8];
            app.ImportrefButton.Position = [24 11 100 22];
            app.ImportrefButton.Text = 'Import ref';

            % Create TissueSaturationTab
            app.TissueSaturationTab = uitab(app.TabGroup);
            app.TissueSaturationTab.Title = 'Tissue Saturation';

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.TissueSaturationTab);
            title(app.UIAxes3, '')
            xlabel(app.UIAxes3, 'Time (s)')
            ylabel(app.UIAxes3, 'StO2 (%)')
            app.UIAxes3.PlotBoxAspectRatio = [1.9453125 1 1];
            app.UIAxes3.Position = [1 1 717 348];

            % Create InitialiseButton
            app.InitialiseButton = uibutton(app.TissueSaturationTab, 'push');
            app.InitialiseButton.ButtonPushedFcn = createCallbackFcn(app, @InitialiseButtonPushed, true);
            app.InitialiseButton.BackgroundColor = [0.8 0.8 0.8];
            app.InitialiseButton.Position = [674 205 100 22];
            app.InitialiseButton.Text = 'Initialise';

            % Create nEditFieldLabel
            app.nEditFieldLabel = uilabel(app.TissueSaturationTab);
            app.nEditFieldLabel.HorizontalAlignment = 'right';
            app.nEditFieldLabel.Position = [674 164 49 22];
            app.nEditFieldLabel.Text = 'n';

            % Create nEditField
            app.nEditField = uieditfield(app.TissueSaturationTab, 'numeric');
            app.nEditField.Limits = [1 1000];
            app.nEditField.Position = [736 164 38 22];
            app.nEditField.Value = 10;

            % Create ControlsPanel
            app.ControlsPanel = uipanel(app.UIFigure);
            app.ControlsPanel.Title = 'Controls';
            app.ControlsPanel.Position = [6 53 136 372];

            % Create StartButton
            app.StartButton = uibutton(app.ControlsPanel, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.BackgroundColor = [0.4706 0.6706 0.1882];
            app.StartButton.Position = [18 205 100 22];
            app.StartButton.Text = 'Start';

            % Create MarkEventButton
            app.MarkEventButton = uibutton(app.ControlsPanel, 'push');
            app.MarkEventButton.ButtonPushedFcn = createCallbackFcn(app, @MarkEventButtonPushed, true);
            app.MarkEventButton.BackgroundColor = [0.9294 0.6902 0.1294];
            app.MarkEventButton.Position = [17 90 100 22];
            app.MarkEventButton.Text = 'Mark Event';

            % Create StopButton
            app.StopButton = uibutton(app.ControlsPanel, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.BackgroundColor = [0.851 0.3294 0.102];
            app.StopButton.Position = [18 169 100 22];
            app.StopButton.Text = 'Stop';

            % Create EventDetailsEditFieldLabel
            app.EventDetailsEditFieldLabel = uilabel(app.ControlsPanel);
            app.EventDetailsEditFieldLabel.HorizontalAlignment = 'right';
            app.EventDetailsEditFieldLabel.Position = [25 66 76 22];
            app.EventDetailsEditFieldLabel.Text = 'Event Details';

            % Create EventDetailsEditField
            app.EventDetailsEditField = uieditfield(app.ControlsPanel, 'text');
            app.EventDetailsEditField.Position = [17 37 100 22];

            % Create LogeventdetailsButton
            app.LogeventdetailsButton = uibutton(app.ControlsPanel, 'push');
            app.LogeventdetailsButton.ButtonPushedFcn = createCallbackFcn(app, @LogeventdetailsButtonPushed, true);
            app.LogeventdetailsButton.BackgroundColor = [0.9294 0.6902 0.1294];
            app.LogeventdetailsButton.Position = [15 11 106 22];
            app.LogeventdetailsButton.Text = 'Log event details';

            % Create miniCYRILatColumbiaLabel
            app.miniCYRILatColumbiaLabel = uilabel(app.UIFigure);
            app.miniCYRILatColumbiaLabel.FontName = 'Lucida Sans Typewriter';
            app.miniCYRILatColumbiaLabel.FontSize = 18;
            app.miniCYRILatColumbiaLabel.FontWeight = 'bold';
            app.miniCYRILatColumbiaLabel.Position = [368 434 275 33];
            app.miniCYRILatColumbiaLabel.Text = 'miniCYRIL at Columbia';

            % Create LiveTemperatureCEditFieldLabel
            app.LiveTemperatureCEditFieldLabel = uilabel(app.UIFigure);
            app.LiveTemperatureCEditFieldLabel.HorizontalAlignment = 'right';
            app.LiveTemperatureCEditFieldLabel.Position = [719 14 118 22];
            app.LiveTemperatureCEditFieldLabel.Text = 'Live Temperature (C)';

            % Create LiveTemperatureCEditField
            app.LiveTemperatureCEditField = uieditfield(app.UIFigure, 'numeric');
            app.LiveTemperatureCEditField.Limits = [-40 40];
            app.LiveTemperatureCEditField.Position = [852 14 44 22];
            app.LiveTemperatureCEditField.Value = -15;

            % Create LiveTimesEditFieldLabel
            app.LiveTimesEditFieldLabel = uilabel(app.UIFigure);
            app.LiveTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.LiveTimesEditFieldLabel.Position = [14 14 75 22];
            app.LiveTimesEditFieldLabel.Text = 'Live Time (s)';

            % Create LiveTimesEditField
            app.LiveTimesEditField = uieditfield(app.UIFigure, 'numeric');
            app.LiveTimesEditField.Position = [100 14 42 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = miniCYRIL_at_Columbia_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @getStarted)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end