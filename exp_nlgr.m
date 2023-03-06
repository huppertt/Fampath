%        Name   : Name of the parameter (a string).
%        Unit   : unit of the parameter (a string).
%        Value  : value of the parameter (a finite real scalar, vector or
%                 2-dimensional matrix).
%        Minimum: minimum values of the parameter (a real scalar, a vector or
%                 a 2-dimensional matrix).
%        Maximum: maximum values of the parameter (a real scalar, a vector or
%                 a 2-dimensional matrix).
%        Fixed  : a boolean, a boolean vector or a boolean 2-dimensional
%                 matrix specifying whether the parameter is fixed or not.

Parameters(1).Name='MIE_A';
Parameters(1).Unit='none';
Parameters(1).Value=1.7;
Parameters(1).Minimum=.5;
Parameters(1).Maximum=5;
Parameters(1).Fixed=false;

Parameters(2).Name='MIE_B';
Parameters(2).Unit='none / 1000';
Parameters(2).Value=24;
Parameters(2).Minimum=0;
Parameters(2).Maximum=50;
Parameters(2).Fixed=false;

Parameters(3).Name='WF';
Parameters(3).Unit='%';
Parameters(3).Value=48;
Parameters(3).Minimum=30;
Parameters(3).Maximum=54;
Parameters(3).Fixed=false;

Parameters(4).Name='ox-redCCO (OD/M/cm)_Moody';
Parameters(4).Unit='uM';
Parameters(4).Value=.8;
Parameters(4).Minimum=0;
Parameters(4).Maximum=4;
Parameters(4).Fixed=false;


Parameters(5).Name='oxCCO (OD/M/cm)_Moody';
Parameters(5).Unit='uM';
Parameters(5).Value=.8;
Parameters(5).Minimum=0;
Parameters(5).Maximum=4;
Parameters(5).Fixed=false;


Parameters(6).Name='redCCO (OD/M/cm)_Moody';
Parameters(6).Unit='uM';
Parameters(6).Value=.8;
Parameters(6).Minimum=0;
Parameters(6).Maximum=4;
Parameters(6).Fixed=false;


Parameters(7).Name='Fat Soybean (OD/cm)_Ekker';
Parameters(7).Unit='%';
Parameters(7).Value=.1;
Parameters(7).Minimum=0;
Parameters(7).Maximum=4;
Parameters(7).Fixed=false;

Parameters(8).Name='DC';
Parameters(8).Unit='%';
Parameters(8).Value=0;
Parameters(8).Minimum=-3E3;
Parameters(8).Maximum=1E4;
Parameters(8).Fixed=false;

Parameters(9).Name='gain_hbo';
Parameters(9).Unit='uM';
Parameters(9).Value=0;
Parameters(9).Minimum=-50;
Parameters(9).Maximum=50;
Parameters(9).Fixed=false;

Parameters(10).Name='gain_hb';
Parameters(10).Unit='uM';
Parameters(10).Value=0;
Parameters(10).Minimum=-50;
Parameters(10).Maximum=50;
Parameters(10).Fixed=false;

Parameters(11).Name='gain_cco_redox';
Parameters(11).Unit='uM';
Parameters(11).Value=0;
Parameters(11).Minimum=-50;
Parameters(11).Maximum=50;
Parameters(11).Fixed=false;



InitialStates(1).Name='HbO';
InitialStates(1).Unit='uM';
InitialStates(1).Value=50;
InitialStates(1).Minimum=10;
InitialStates(1).Maximum=75;
InitialStates(1).Fixed=false;

InitialStates(2).Name='Hb';
InitialStates(2).Unit='uM';
InitialStates(2).Value=20;
InitialStates(2).Minimum=10;
InitialStates(2).Maximum=35;
InitialStates(2).Fixed=false;

InitialStates(3).Name='CCO_red';
InitialStates(3).Unit='uM';
InitialStates(3).Value=1;
InitialStates(3).Minimum=.1;
InitialStates(3).Maximum=10;
InitialStates(3).Fixed=false;

lstValidLambda=find(data.Wavelengths>700 & data.Wavelengths<950);


Ext=tissue_specific_extinction_coefficient_650to1000;

% ---------------------------------------------------------------------------------------------------------
% 1-Wavelengths(nm)	
% 2-Water(OD/M/cm)	
% 3-HbO2 (OD/M/cm)
% 4-HHb (OD/M/cm)
% 5-Difference Cytochrome Oxidase (OD/M/cm) 
% 6-ox-redCCO (OD/M/cm)_Moody 
% 7-oxCCO (OD/M/cm)_Moody 
% 8-redCCO (OD/M/cm)_Moody 
% 9-Fat Soybean (OD/cm)_Ekker
% ---------------------------------------------------------------------------------------------------------

% lambda, HbO2, Hb, diffCCO, WF, CCO,
Ext=Ext(:,[1 3 4 5 2 6 7 8 9]);

Ext_int=zeros(length(lstValidLambda),size(Ext,2));
Ext_int(:,1)=data.Wavelengths(lstValidLambda);
for i=2:size(Ext,2)
    Ext_int(:,i)=interp1(Ext(:,1),Ext(:,i),Ext_int(:,1),'linear','extrap');
end

AdditionalData.Lambda=data.Wavelengths(lstValidLambda)'*1E-9;
AdditionalData.ext_coef=Ext_int(:,2:end);
AdditionalData.rho=3.2;
AdditionalData.Reff=0.493;


refd=ted_BH.raw_data;
refd=nirs.math.tddr(refd',5)';
refd=nirs.math.tddr(refd,5);
%I = A*Io*exp(-mua)
% mua = log(A)+log(I0)-log(I)
refd = log(ones(size(refd,1),1)*abs(ref.light)')-log(refd)+log(10);

mrefd=median(refd,1);
refd=(refd'-PCA*inv(PCA'*PCA)*PCA'*refd')';
reflect=refd+ones(size(refd,1),1)*mrefd;
% 
% reflect = (data.Spectra(:,lstValidLambda)./(ones(size(data.Spectra,1),1)*data.Ref(lstValidLambda)));
% u=data.Events(:,1);
z = iddata(reflect(:,lstValidLambda), u, median(diff(data.Time)));

Order         = [length(lstValidLambda) 1 3];  % Model orders [ny nu nx].
nlgr = idnlgrey('hyperspectral_m', Order, Parameters, InitialStates, 0, ...
    'Name', 'Hyperspectral fit');

nlgr.FileArgument={AdditionalData};


% C. Estimate the parameters of the DC-motor model object.
opt = nlgreyestOptions('Display','on');
nlgr = nlgreyest(z, nlgr, opt);










for i=1:length(paramorder)
    id=find(ismember({obj.parameters.Name},paramorder{i}));
    pValues(i)=obj.InitialStates(id).value;
    pValuesMin{i}=obj.InitialStates(id).range(1);
    pValuesMax{i}=obj.InitialStates(id).range(2);
    pValuesFixed{i}=obj.InitialStates(id).fixed;
end

FileName      = @nirs.vascular.models.WKM.WKM_dyn;   % File describing the WKM model structure.
Order         = [2 2 7];   % Model orders [ny nu nx].
Parameters    = pValues;   % Initial parameters. Np = 7.
InitialStates = [0; 0;1;1;1;1;1];   % Initial initial states.
Ts            = 0;           % Time-continuous system.
nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
    'Name', 'Windkessel model');

set(nlgr, 'InputName', {'flow-inducing','CMRO2-inducing'}, 'InputUnit', {'%','%'},...
    'OutputName', {'HbO2', 'HbR'}, ...
    'OutputUnit', {'uM', 'uM'},                         ...
    'TimeUnit', 's');

nlgr = setinit(nlgr,'Name', obj.states);
nlgr = setpar(nlgr, 'Name', paramorder);
nlgr = setpar(nlgr, 'Minimum', pValuesMin);
nlgr = setpar(nlgr, 'Maximum', pValuesMax);
nlgr = setpar(nlgr, 'Fixed', pValuesFixed);

nlgr.SimulationOptions.AbsTol = 1e-6;
nlgr.SimulationOptions.RelTol = 1e-5;