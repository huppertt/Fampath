function FAMPATH_estimate_model(data)

% data.Spectra
% data.Reference
% data.Wavelengths
% data.Time
% data.Events
% 
% Events=Dictionary;
% 
% stim=nirs.design.StimulusEvents;
% stim.name='BreathHold';
% stim.onset=[0 60 120 180 240]'+tmp.Time(min(find(tmp.Events(:,1)>.5)));
% stim.dur=[30 30 30 30 30]';
% stim.amp=[1 1 1 1 1]';
% Events('BreathHold')=stim;
% 
% stim.name='Rapid';
% stim.onset=[0 60 120 180 240]'+tmp.Time(min(find(tmp.Events(:,1)>.5)))+570;
% stim.dur=[30 30 30 30 30]';
% stim.amp=[1 1 1 1 1]';
% Events('Rapid')=stim;


options.wavelength_range=[650 1000];
options.smoother=@(x)sgolayfilt(x,4,9,[],2);
options.Reff=0.493;
options.rho=37;

validwavelengths=find(data.Wavelength>=options.wavelength_range(1) &...
    data.Wavelength<=options.wavelength_range(2));

options.wavelengths=data.Wavelength(validwavelengths);

c=tissue_specific_extinction_coefficient_650to1000;
for i=2:size(c,2) 
    coef(:,i-1)=interp1(c(:,1),c(:,i),data.Wavelength(validwavelengths),'linear','extrap'); 
end;
options.coeff=coef;

ntps=length(data.Time);
nmeas=length(validwavelengths);

Y=data.Spectrum';
Y=medfilt2(Y,[5 5]);
Y=Y-ones(size(Y,1),1)*median(Y(1:20,:),1);


Y=Y(validwavelengths,:);
Ref=interp1(data.calib_Wavelength,data.Lamp,options.wavelengths)';
DarkNoise=interp1(data.calib_Wavelength,data.DarkNoise,options.wavelengths)';




meas.reflectance = options.smoother(Y'./(ones(ntps,1)*options.smoother(Ref')));
meas.diff_1_reflectance = options.smoother(diff(meas.reflectance,1,2));
meas.diff_2_reflectance = options.smoother(diff(meas.diff_1_reflectance,1,2));

meas.refl_whitening=diag(1./sqrt(var(log(meas.reflectance),[],1)));
meas.diff1_whitening=diag(1./sqrt(var(log(meas.diff_1_reflectance),[],1)));
meas.diff2_whitening=diag(1./sqrt(var(log(meas.diff_2_reflectance),[],1)));

% 
% Boundaries = [95, 80, 60, 30, 5;...
%               20, 0, 0, .05,0];

fitoptions = optimoptions('lsqnonlin','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-16,'MaxFunctionEvaluations',1000,'MaxIterations',1000,'Display','iter',...
    'Algorithm','levenberg-marquardt');

% 
% pp=mean(Boundaries,1);
% rBoundaries=Boundaries;

y1=median(meas.reflectance,1);
y2=median(meas.diff_1_reflectance,1);
y3=median(meas.diff_2_reflectance,1);
y4=[y1 y2 y3];

costfun1=@(param)(abs(log(sim_refl_model(param,options.wavelengths,options,'refl'))-log(y1))');
costfun2=@(param)(abs(log(sim_refl_model(param,options.wavelengths,options,'diff_1'))-log(y2))');
costfun3=@(param)(abs(log(sim_refl_model(param,options.wavelengths,options,'diff_2'))-log(y3))');
costfun4=@(param)(abs(log(sim_refl_model(param,options.wavelengths,options,'all'))-log(y4))');


pp=[0 1 1 (pinv(coef'*coef)*coef'*log(median(meas.reflectance,1)'))'];
x=pp;

while(1)
x=lsqnonlin(costfun4,x,[],[],fitoptions);

figure(1);
subplot(4,1,1);
cla;  hold on;
plot(options.wavelengths,costfun1(pp),'k');
plot(options.wavelengths,y1,'r');

subplot(4,1,2);
cla;  hold on;
plot(options.wavelengths(1:end-1),costfun2(pp),'k');
plot(options.wavelengths(1:end-1),y2,'r');

subplot(4,1,3);
cla;  hold on;
plot(options.wavelengths(1:end-2),costfun3(pp),'k');
plot(options.wavelengths(1:end-2),y3,'r');

subplot(4,1,4);
cla;  hold on;
plot(costfun4(pp),'k');
plot(y4,'r');
pause;
end



