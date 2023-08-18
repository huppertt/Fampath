classdef Physiol < handle
    properties
        pulseox_connected=false;
        capno_connected=false;
        capno_port='COM3';
        
    end
    properties(Hidden=true)
        ble_channel;
        ble_device;
        capno_device;
        
        data_buffer=1000;
        
        raw_pulse_data;
        raw_co2_data;
        raw_EtCO2;
        raw_RespRate;
        raw_PulseRate;
        raw_SpO2;
        
        raw_pulse_time;
        raw_SpO2_time;
        raw_EtCO2_time;
        raw_CO2_time;
        
        raw_pulse_cnt;
        raw_SpO2_cnt;
        raw_EtCO2_cnt;
        raw_CO2_cnt;
        isrunning=false;
        start_time;
        last_sample;
    end
    
    methods
        function obj=Physiol()
            %initializer
            
            % find the BLE pulse oximeter
            try
                ble_devicelist=blelist;
                chIdx=find(contains(ble_devicelist.Name,'OxySmart'));
                if(isempty(chIdx))
                    warning('Unable to find OxySmart PulseOx');
                    obj.pulseox_connected=false;
                else
                    try
                        
                        obj.ble_device=ble(ble_devicelist.Address(chIdx));
                        obj.pulseox_connected=true;
                        i=find((contains(obj.ble_device.Characteristics.ServiceUUID,'-')), 1, 'last' );
                        obj.ble_channel=characteristic(obj.ble_device,...
                            obj.ble_device.Characteristics.ServiceUUID(i),...
                            obj.ble_device.Characteristics.CharacteristicUUID(i));
                        obj.ble_channel.DataAvailableFcn=@(d,t)local_add_pulse_data(obj,d,t);
                    catch
                        warning('Unable to connect with OxySmart PulseOx');
                        obj.pulseox_connected=false;
                    end
                end
            end
            try
                obj.capno_device=serial(obj.capno_port);
                set(obj.capno_device,'BaudRate',9600,'DataBits',8,...
                    'FlowControl','none','Parity','none');
                set(obj.capno_device,'InputBufferSize',obj.data_buffer*8);
                set(obj.capno_device,'TimerFcn',@(obj,event,varargin)add_capno_data_local(obj));
                set(obj.capno_device,'TimerPeriod',.25);
                set(obj.capno_device,'UserData',obj);
                obj.capno_connected=true;
            catch
                obj.capno_connected=false;
            end
            
            
            obj.reset();
            obj.isrunning=false;
            
            
            
        end
        
        
        function local_add_pulse_data(obj,device,events)
            if(~obj.isrunning)
                return;
            end
            try
                [d,t]=read(device,'oldest');
                t=datenum(t);
                if(d(3)==15 & d(4)==8)
                    obj.raw_SpO2_cnt=obj.raw_SpO2_cnt+1;
                    obj.raw_PulseRate(obj.raw_SpO2_cnt)=d(7);
                    obj.raw_SpO2(obj.raw_SpO2_cnt)=d(6);
                    obj.raw_SpO2_time(obj.raw_SpO2_cnt)=t;
                elseif(d(3)==15 & d(4)==7)
                    obj.raw_pulse_cnt=obj.raw_pulse_cnt+1;
                    b=(d(11)+d(10)*2^8+d(9)*2^16+...
                        d(8)*2^24+d(7)*2^32+d(8)*2^40)/(2^48);
                    obj.raw_pulse_data(obj.raw_pulse_cnt)=b;
                    obj.raw_pulse_time(obj.raw_pulse_cnt)=t;
                end
            end
        end
            
        
        function start(obj)
            obj.isrunning=true;
            fopen(obj.capno_device);
            obj.start_time=now;
            obj.last_sample=0;
        end
        
        function stop(obj)
            obj.isrunning=false;
            fclose(obj.capno_device);
        end
        
        function data = getdata(obj)
            data=struct;
            data.pulse=obj.raw_pulse_data(1:obj.raw_pulse_cnt);
            data.time=obj.raw_pulse_time(1:obj.raw_pulse_cnt);
            data.PulseRate=obj.raw_PulseRate(1:obj.raw_SpO2_cnt);
            data.SpO2=obj.raw_SpO2(1:obj.raw_SpO2_cnt);
            data.SpO2_time=obj.raw_SpO2_time(1:obj.raw_SpO2_cnt);
            
            data.EtCO2_time=obj.raw_EtCO2_time(1:obj.raw_EtCO2_cnt);
            data.EtCO2=obj.raw_EtCO2(1:obj.raw_EtCO2_cnt);
            data.RespRate=obj.raw_RespRate(1:obj.raw_EtCO2_cnt);
            
            data.co2=obj.raw_co2_data(1:obj.raw_CO2_cnt);
            
            reset(obj);
        end
        
        
        
        
        function reset(obj)
            obj.raw_pulse_data=zeros(obj.data_buffer,1);
            obj.raw_co2_data=zeros(obj.data_buffer,1);
            obj.raw_EtCO2=zeros(obj.data_buffer,1);
            obj.raw_RespRate=zeros(obj.data_buffer,1);
            obj.raw_PulseRate=zeros(obj.data_buffer,1);
            obj.raw_SpO2=zeros(obj.data_buffer,1);
            obj.raw_pulse_time=zeros(obj.data_buffer,1);
            obj.raw_SpO2_time=zeros(obj.data_buffer,1);
            obj.raw_CO2_time=zeros(obj.data_buffer,1);
            obj.raw_EtCO2_time=zeros(obj.data_buffer,1);
            obj.raw_pulse_cnt=0;
            obj.raw_SpO2_cnt=0;
            obj.raw_CO2_cnt=0;
            obj.raw_EtCO2_cnt=0;
        end
        
        
        function add_capno_data(obj,varargin)
            
            try
                data=fread(obj.capno_device,obj.capno_device.BytesAvailable);
                time_elapsed=now-obj.start_time;
                
                EtCO2=[];
                CO2=[];
                RR=[];
                
                startIndices=find(data==hex2dec('FF'));
                for i=1:length(startIndices)
                    try
                        local_data=data(startIndices(i)+[0:5]);
                        
                        keybyte=local_data(2);
                        if(keybyte==170)
                            % 0xAA -CO2 in mmHg
                            scale=1;
                        else
                            %0xAB - CO2 in kPA
                            scale=1/10;
                        end
                        
                        
                        IDbyte=local_data(3);
                        switch(IDbyte)
                            case 240
                                %0xF0- EtCO2.  2 times per sec
                                EtCO2(end+1)=scale*local_data(4);
                            case 241
                                %0xF1- Resp Rate.  2 times per sec
                                RR(end+1)=local_data(4);
                            case 242
                                %0xF2- Co2 graph.  4 times per sec
                                CO2(end+1)=scale*local_data(4);
                            case 243
                                %0xF3- capno alarm
                            case 248
                                %0xF8- system alarm
                            case 250
                                %0xFA- limit alarm
                            case 253
                                %0xFD- alarm status
                        end
                    end
                    %     checksum=256*local_data(5)+local_data(6);
                    %     checksum2=local_data(3)+local_data(4);
                    %     if(checksum~=checksum2)
                    %         error('checksum missmatch');
                    %     end
                end
                
                n=min([length(EtCO2) length(RR)]);
                EtCO2=EtCO2(1:n);
                RR=RR(1:n);
                
                cnt=obj.raw_EtCO2_cnt+1;
                obj.raw_EtCO2_cnt=cnt+length(EtCO2)-1;
                
                t=linspace(obj.last_sample,time_elapsed,length(EtCO2)+1);
                obj.last_sample=time_elapsed;
                obj.raw_EtCO2_time(cnt:obj.raw_EtCO2_cnt)=t(2:end)*86408;
                obj.raw_EtCO2(cnt:obj.raw_EtCO2_cnt)=EtCO2;
                obj.raw_RespRate(cnt:obj.raw_EtCO2_cnt)=RR;
                
                
                cnt=obj.raw_CO2_cnt+1;
                obj.raw_CO2_cnt=cnt+length(CO2)-1;
                obj.raw_co2_data(cnt:obj.raw_CO2_cnt)=CO2;
                
            end
            
        end
    end
end


function add_capno_data_local(obj,varargin)
     obj=obj.UserData;
     try
        obj.add_capno_data;
     end
end