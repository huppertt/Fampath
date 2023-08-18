function varargout=mini_CYRIL_update_data(varargin)
app=get(findall(0,'Type','figure','Name','MINICYRIL'),'userdata');

% Step 1- grab new data and update the spectrum display
[data,wavelength]=app.system_info.Spectrometer_driver.get_spectrum();
app.LiveTimesEditField.Value= 24*3600*(now-app.daq_start_time);

if(~isempty(app.MarkEventButton.UserData))
    fprintf(app.file_info.file_fid,'Event %s\n',app.MarkEventButton.UserData);
    app.MarkEventButton.UserData=[];
    app.stored_data.events(app.stored_data.count)=999;
    yaxis=get(app.UIAxes_Concentration,'Ylim');
    t=24*3600*(now-app.daq_start_time);
    plot(app.UIAxes_Concentration,[t t],[yaxis]);
    %TODO- add event line on plot
 end
% 
% if(~isempty(app.system_info.inlet))
%     try
%         [mrks,ts] = app.system_info.inlet.pull_sample(.05);
%         if(~isempty(ts))
%             app.stored_data.events(app.stored_data.count)=mrks;
%             yaxis=get(app.UIAxes_Concentration,'Ylim');
%             t=24*3600*(now-app.daq_start_time);
%             plot(app.UIAxes_Concentration,[t t],[yaxis]);
%         end
%     end
% end

maxtime=app.DisplayWindowsEditField.Value;


app.stored_data.time(app.stored_data.count)=app.LiveTimesEditField.Value;
app.stored_data.raw_data(app.stored_data.count,:)=data;
app.stored_data.count=app.stored_data.count+1;


set(app.drawing_handles.spectrum.spectrum,'Ydata',data,'Xdata',wavelength);
xlim(app.UIAxes_Spectra, [600 950])
drawnow;
time=datestr(now);

if(app.saving_data)
    fprintf(app.file_info.file_fid,'%s,',time);
    fprintf(app.file_info.file_fid,'%d,',data);
    fprintf(app.file_info.file_fid,'\n');
end
try


app.conc_changes.update(data);

HbO2=app.conc_changes.HBO2;
HbR=app.conc_changes.HHB;
CCO=app.conc_changes.CCO;
app.stored_data.conc_data(app.stored_data.count,:)=[HbO2 HbR CCO];

% Update the conc time-traces
xdata=get(app.drawing_handles.conc.hbo,'Xdata');
xdata=[xdata(:)' app.LiveTimesEditField.Value];
ydata=get(app.drawing_handles.conc.hbo,'Ydata');
set(app.drawing_handles.conc.hbo,'Ydata',[ydata(:)' HbO2],'Xdata',xdata);
ydata=get(app.drawing_handles.conc.hbr,'Ydata');
set(app.drawing_handles.conc.hbr,'Ydata',[ydata(:)' HbR],'Xdata',xdata);
ydata=get(app.drawing_handles.conc.cco,'Ydata');
set(app.drawing_handles.conc.cco,'Ydata',[ydata(:)' CCO],'Xdata',xdata);


 
% Finally update the Capnometry and Pulse ox traces
try
    data=app.physiol_monitor.getdata();
    app.stored_data.phys_data{app.stored_data.count}=data;
    if(~isempty(data.EtCO2))
        app.Label_ETCO2_Live.Text=num2str(data.EtCO2(end));
        app.Label_RR_Live.Text=num2str(data.RespRate(end));

        xdata=get(app.drawing_handles.physio.etco2,'Xdata');
        xdata=[xdata(:)' app.LiveTimesEditField.Value];
        ydata=get(app.drawing_handles.physio.etco2,'Ydata');
        set(app.drawing_handles.physio.etco2,'Ydata',[ydata(:)' data.EtCO2(end)],'Xdata',xdata);
    end
    if(~isempty(data.pulse))
        xdata=get(app.drawing_handles.physio.SpO2,'Xdata');
        if(length(xdata)==1); 
            set(app.drawing_handles.physio.SpO2,'Userdata',data.time(1));
        end
        data.time=(data.time-get(app.drawing_handles.physio.SpO2,'Userdata'));
        xdata=[xdata(:)' data.time(:)'];
        ydata=get(app.drawing_handles.physio.SpO2,'Ydata');
        ydata(ydata==0)=NaN;
        set(app.drawing_handles.physio.SpO2,'Ydata',[ydata(:)' data.pulse(:)'],'Xdata',xdata);
        
    end
    
    if(~isempty(data.PulseRate))
        app.Label_Pulse_Live.Text=num2str(data.PulseRate(end));
        app.Label_SaO2_Live.Text=num2str(data.SpO2(end));
    end
end


tt=app.LiveTimesEditField.Value;
xl(1)=max(tt-maxtime,0);
xl(2)=tt;
app.UIAxes_Concentration.XLim=xl;
app.UIAxes_EtCO2.XLim=xl;
app.UIAxes_SO2.XLim=xl;
drawnow;

%app.LiveTemperatureCEditField.Value = double(spectrometer.detectorTemperatureDegC); %display temp

refreshdata(app.MINICYRILUIFigure);
catch
    disp(lasterr);
end
