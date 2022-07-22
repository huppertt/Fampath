classdef Spectrometer < handle
    % Class to control the FAMPATH spectrometer
    
    % This example requires the following:
    % * A 64-bit Microsoft(R) Windows(R)
    % * Ocean Optics spectrometer USB2000
    % * Install OmniDriver downloadable from http://www.oceanoptics.com/
    % * OmniDriver.mdd available from MATLAB Central
    
    
    properties       
        IntegrationTime =1000;
        BoxcarWidth=40;
        ScansToAverage=1;
        CorrectForDetectorNonlinearity=false
        CorrectForElectricalDark=false;
        detectorTECSetpointDegC;
        system_info=[];

    end
    
    properties(Dependent=true)
        model_name;
        model_serial;
        TEC_temperature;
    end
    
    properties(Hidden=true)
        spectrometerObj
        % Spectrometer index to use (first spectrometer by default).
        spectrometerIndex = 0;
        % Channel index to use (first channel by default).
        channelIndex = 0;
        is_OceanOptics=true;
    end
    
    methods
        
     
        function obj=Spectrometer(is_OceanOptics)
            obj.is_OceanOptics=is_OceanOptics;
            obj.spectrometerIndex=0;
        end
        
        function temp = get.TEC_temperature(obj)
         if(~obj.is_OceanOptics)
             temp=obj.spectrometer.detectorTemperatureDegC
         else
             temp=NaN;
         end
        end
        
        function name = get.model_name(obj)
            name = obj.system_info(obj.spectrometerIndex+1).spectrometerName;
        end
        
        function serial = get.model_serial(obj)
            serial = obj.system_info(obj.spectrometerIndex+1).spectrometerSerialNumber;
        end
        
        
        function obj = set.IntegrationTime(obj,IntegrationTime)
            obj.IntegrationTime=IntegrationTime;
            if(obj.is_OceanOptics)
                invoke(obj.spectrometerObj, 'setIntegrationTime', obj.spectrometerIndex, 0, 1000*obj.IntegrationTime/obj.ScansToAverage);
          else
                 obj.spectrometerObj.integrationTimeMSobj.IntegrationTime;
      
            end
        end
        
        function obj = set.CorrectForDetectorNonlinearity(obj,enable)
            obj.CorrectForDetectorNonlinearity=enable;
            if(obj.is_OceanOptics)
                invoke(obj.spectrometerObj, 'setCorrectForDetectorNonlinearity',...
                    obj.spectrometerIndex, obj.channelIndex,obj.CorrectForDetectorNonlinearity);
            end
        end
        
        function obj = set.CorrectForElectricalDark(obj,enable)
            obj.CorrectForElectricalDark=enable;
            if(obj.is_OceanOptics)
                invoke(obj.spectrometerObj, 'setCorrectForElectricalDark',...
                    obj.spectrometerIndex, obj.channelIndex, obj.CorrectForElectricalDark);
            end
        end
        
        function set.detectorTECSetpointDegC(obj,detectorTECSetpointDegC)
            if(~obj.is_OceanOptics)
                obj.detectorTECSetpointDegC=detectorTECSetpointDegC;
                obj.spectrometer.detectorTECEnabled = 1;
                obj.spectrometer.detectorTECSetpointDegC = detectorTECSetpointDegC;
            end
        end
        
        function set.spectrometerIndex(obj,spectrometerIndex)
            obj.spectrometerIndex=spectrometerIndex;
            if(~obj.is_OceanOptics)
                obj.spectrometerObj = DLL.getSpectrometer(obj.spectrometerIndex);
            else
            end
        end
        
        function [data,wavelengths]=get_spectrum(obj)
            %% Aquire the spectrum.
            if(obj.is_OceanOptics)
            wavelengths = invoke(obj.spectrometerObj,...
                'getWavelengths', obj.spectrometerIndex, obj.channelIndex);
            % Get the wavelengths of the first spectrometer and save them in a double
            % array.
            
                data = invoke(obj.spectrometerObj, 'getSpectrum', obj.spectrometerIndex);
            else
                wavelengths = obj.spectrometer.wavelengths;
                data = obj.spectrometer.getSpectrum();
            end
        end
        
        
        function initialize(obj)
            if(obj.is_OceanOptics)
                obj.spectrometerObj = icdevice('OceanOptics_OmniDriver.mdd');
                connect(obj.spectrometerObj);
                
                
                % Get number of spectrometers connected.
                numOfSpectrometers = invoke(obj.spectrometerObj, 'getNumberOfSpectrometersFound');
                
                for idx=1:numOfSpectrometers
                    info.spectrometerName = invoke(obj.spectrometerObj, 'getName', idx-1);
                    info.spectrometerSerialNumber = invoke(obj.spectrometerObj, 'getSerialNumber', idx-1);
                    if(isempty(obj.system_info))
                        obj.system_info=info;
                    else
                        obj.system_info(idx,1)=info;
                    end
                end
                invoke(obj.spectrometerObj, 'setIntegrationTime', obj.spectrometerIndex, 0, 1000*obj.IntegrationTime/obj.ScansToAverage);
                invoke(obj.spectrometerObj,'setBoxcarWidth',obj.spectrometerIndex,0,obj.BoxcarWidth) 
                invoke(obj.spectrometerObj, 'setCorrectForElectricalDark',...
                    obj.spectrometerIndex, obj.channelIndex, obj.CorrectForElectricalDark);
                invoke(obj.spectrometerObj, 'setCorrectForDetectorNonlinearity',...
                    obj.spectrometerIndex, obj.channelIndex,obj.CorrectForDetectorNonlinearity);
                invoke(obj.spectrometerObj,'setScansToAverage',obj.spectrometerIndex,obj.ScansToAverage) 
                 
            else
                dll = NET.addAssembly('C:\Program Files\Wasatch Photonics\Wasatch.NET\WasatchNET.dll');
                DLL = WasatchNET.Driver.getInstance();
                %open spectrometer and set temperature
                numOfSpectrometers = DLL.openAllSpectrometers();
                
                for idx=1:numOfSpectrometers
                    obj.spectrometerObj = DLL.getSpectrometer(idx-1);
                    info.spectrometerName = char(spectrometer.model);
                    info.spectrometerSerialNumber =  char(spectrometer.serialNumber);
                    if(isempty(obj.system_info))
                        obj.system_info=info;
                    else
                        obj.system_info(idx,1)=info;
                    end
                end
                obj.spectrometerObj = DLL.getSpectrometer(obj.spectrometerIndex);
         

            end
            
        end
        
        function delete(obj)
            if(obj.is_OceanOptics)
                try
                    disconnect(obj.spectrometerObj);
                end
                try
                    delete(obj.spectrometerObj);
                end
            end
        end
        
        
    end
end

