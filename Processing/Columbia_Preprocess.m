function [NewData2,Absf]=Columbia_Preprocess(subjid,timingfile)

PL_method=2; %1=DPF, 2=Water fitting

if(nargin<2)
    timingfile=[];
end

ProcessingNotes={};
ProcessingNotes{end+1}=['Analysis run on ' datestr(now)];
waverange=[645 917];

filemat=rdir(['**/*' subjid '*.mat']);
filexls=dir(['**/*' subjid '*Capno.xls*']);

if(length(filemat)==0)
    warning('unable to find MAT file');
    return
end
data=load(filemat(1).name);
ProcessingNotes{end+1}=['Read data from file:'  filemat(1).name];


if(isfield(data,'Spectra'))
    ProcessingNotes{end+1}='Old FAMPATH data file detected';
    disp(['older data version detected']);
    lstW=find(data.Wavelengths>=waverange(1) & data.Wavelengths<=waverange(2));
    OD = data.Spectra(:,lstW);
    time=data.Time;
    if(length(filexls)==0)
        warning('unable to find CAPNO XLS file');
        return
    end
    tbl=readtable(filexls(1).name,'Sheet','Data');
    oldver=true;
else
    ProcessingNotes{end+1}='New FAMPATH data file detected';
    lstW=find(data.wavelengths>=waverange(1) & data.wavelengths<=waverange(2));

    OD = data.raw_data(:,lstW);
    time=data.time;
    tbl=data.Phys;
    oldver=false;
end

warning('off','MATLAB:table:ModifiedAndSavedVarnames')

if(~isempty(timingfile))
    ProcessingNotes{end+1}=['Supplimental timing file used: ' timingfile];
    colnames={'ID','Study','Date','Group','Counterbalancing','Session_time',...
        'Experimenter','NIRS_time','Eprime_time','NIRS_offset','Cap_offset','NIRS_errors',...
        'Eprime_errors','Capnometry_errors','Other_errors','Behavior_notes',...
        'Head_circum'};
    timingfile=readtable(timingfile);
    timingfile(:,length(colnames)+1:end)=[];
    timingfile.Properties.VariableNames=colnames;
end

lst=min(find(mean(OD,2)>5*median(mean(OD,2))))-5:length(time);
time(lst)=[];
OD(lst,:)=[];

time_ideal=time(1):time(end);
OD_ideal = zeros(length(time_ideal),size(OD,2));
for i=1:size(OD,2);
    OD_ideal(:,i)=interp1(time,OD(:,i),time_ideal,"linear");
end

[U,S,V]=nirs.math.mysvd(OD_ideal);
U=nirs.math.tddr(U,1,true);
DC=median(U,1);
[fa,fb]=butter(4,[1/300 .4]/2);
U=filtfilt(fa,fb,U);
U=U+ones(size(U,1),1)*DC;

OD_ideal_filtered=U*S*V';

DC=median(OD_ideal,1);
[fa,fb]=butter(4,[1/300 .4]/2);
OD_ideal=filtfilt(fa,fb,OD_ideal)+ones(size(U,1),1)*DC;

ProcessingNotes{end+1}=['Data filtered with 1/300-0.4Hz bandpass filter'];
ProcessingNotes{end+1}=['Motion correction done using PCA-TDDR'];

if(oldver)
    temp_ideal=interp1(data.Time,data.Temp,time_ideal,'linear');

    tbl.Time=tbl.Time*(60*60*24);
    if(length(tbl.Time)~=length(unique(tbl.Time)))
        tbl.Time(:)=tbl.Time(1)+[0:height(tbl)-1]*.25;
    end


    Events=time(data.Events(:,2)>0);
    Events(1)=[];


    if(~isempty(timingfile))
        found=[];
        for i=1:height(timingfile)
            if(contains(subjid,timingfile.ID{i}))
                found=i;
                break;
            end
        end
        if(isempty(timingfile))
            ProcessingNotes{end+1}=['(WARNING) Subject not found in suppl timing file'];
            warning(['unable to find subject ' subjid ' in timing file']);
            timingfile=[];
            info=[];
        else
            info=table2struct(timingfile(found,:));
        end
    end
    % Let's use the fact that EtCO2 becomes invalid during Breath holding to
    % find the sync between the Physiology and the NIRS data
    if(~isempty(info) && ~isnan(info.Cap_offset))
        shift=info.Cap_offset;
        tbl.Time=tbl.Time-shift;
    end


    EtCO2=interp1(tbl.Time,tbl.EtCO2,time_ideal,'linear')';
    RR=interp1(tbl.Time,tbl.RR,time_ideal,'linear')';
    Pulse=interp1(tbl.Time,tbl.Pulse,time_ideal,'linear')';
    SpO2=interp1(tbl.Time,tbl.SpO2,time_ideal,'linear')';

    if(isempty(data.Event_details) & isempty(info))
        ProcessingNotes{end+1}=['(WARNING) Event order not defined'];
        warning('Event details are not labeled')
        Rapid_Times = Events(1)+60*[0:4]';
        Hold_Times  = Events(1+find(diff(Events)>270))+60*[0:4]';
    else

        if(contains(info.Counterbalancing,'R'))
            ProcessingNotes{end+1}=['Using Rapid-Breathing first order'];
            %Start with Rapid
            Rapid_Times = Events(1)+60*[0:4]';
            Hold_Times  = Events(1+find(diff(Events)>270))+60*[0:4]';
        else
            %Start with Hold
            ProcessingNotes{end+1}=['Using Breath-Hold first order'];
            Hold_Times = Events(1)+60*[0:4]';
            Rapid_Times  = Events(1+find(diff(Events)>270))+60*[0:4]';
        end

        wavelengths=data.Wavelengths(lstW);
        referenceOD=data.Ref(lstW);

    end


else
    EtCO2 = interp1(tbl.Capno.time,tbl.Capno.ETCO2,time_ideal,'linear')';
    RR = interp1(tbl.Capno.time,tbl.Capno.RespRate,time_ideal,'linear')';
    if(~isempty(tbl.Pulse.PulseRate))
        Pulse = interp1(tbl.Pulse.time,tbl.Pulse.PulseRate,time_ideal,'linear')';
        SpO2 = interp1(tbl.Pulse.time,tbl.Pulse.SpO2,time_ideal,'linear')';
    else
        ProcessingNotes{end+1}='(WARNING) Pulse oximeter data missing';
        SpO2=nan(size(time_ideal))';
        Pulse=nan(size(time_ideal))';
    end

    if(isempty(find(data.events==1)))
        ProcessingNotes{end+1}=['(WARNING) Eprime marks not found.  Using manual marks'];
        e=time(find(data.events==999));
        Rapid_on=e(1)+60*[0:4]'+300;
        BH_on = Rapid_on(end)+360+60*[0:4]';
    else
        if(length(find(data.events==5))==5)
            a=1;
        else
            a=0;
        end

        Rapid_on=data.time(find(data.events==1+a));
        if(length(Rapid_on)==6)
            Rapid_on(1)=[];
        end
        Rapid_off=data.time(find(data.events==2+a));
        BH_on=data.time(find(data.events==3+a));
        BH_off=data.time(find(data.events==4+a));
    end

    Events(1)=min([Rapid_on; BH_on]);
    Rapid_Times=Rapid_on;
    Hold_Times=BH_on;
    temp_ideal=nan(size(time_ideal));

    wavelengths=data.wavelengths(lstW);
    referenceOD=data.reference.spectrum(lstW);

    info.Study=data.scaninfo.patientID;
    info.Group=data.scaninfo.group;
    info.Date=data.scaninfo.date;
    info.Experimenter='?';
    info.NIRS_errors='NA';
    info.Eprime_errors='NA';
    info.Capnometry_errors='NA';
    info.Other_errors='NA';
    info.Behavior_notes=data.scaninfo.Comments{1};
end



stim=zeros(length(time_ideal),2);
for i=1:length(Hold_Times)
    lst=find(time_ideal>=Hold_Times(i) & time_ideal<=Hold_Times(i)+20);
    stim(lst,1)=i;
end
for i=1:length(Rapid_Times)
    lst=find(time_ideal>=Rapid_Times(i) & time_ideal<=Rapid_Times(i)+20);
    stim(lst,2)=i;
end

if(PL_method==2)
    ProcessingNotes{end+1}='Using Water Fitting method';
else
    ProcessingNotes{end+1}='Using DPF method';
end

[HbO2,Hb,CCO,Abs] = process_Columbia_conc(OD_ideal,wavelengths,referenceOD,PL_method);
[HbO2f,Hbf,CCOf,Absf] = process_Columbia_conc(OD_ideal_filtered,wavelengths,referenceOD,PL_method);

if(PL_method==2)
    [HbO2,Hb,CCO,~] = process_Columbia_conc(OD_ideal,wavelengths,referenceOD,PL_method);
    [HbO2f,Hbf,CCOf,~] = process_Columbia_conc(OD_ideal_filtered,wavelengths,referenceOD,PL_method);
end

%StO2=(HbO2+median(Abs.HbO2))./(Hb+median(Abs.HHb)+HbO2+median(Abs.HbO2));
%StO2f=(HbO2f+median(Absf.HbO2))./(Hbf+median(Absf.HHb)+HbO2f+median(Absf.HbO2));

delete([subjid '_analyzed.xlsx']);

NewData=struct;
NewData.Time=time_ideal'-Events(1)+300;

NewData.Event_Rapid=stim(:,2);
NewData.Event_Hold=stim(:,1);

NewData.HbO2=HbO2;
NewData.Hb=Hb;
NewData.CCO=CCO;
%NewData.StO2=StO2;
NewData.EtCO2=EtCO2;
NewData.RR=RR;
NewData.SpO2=SpO2;
NewData.Pulse=Pulse;
NewData.Temperature=temp_ideal';
NewData=struct2table(NewData);
NewData(NewData.Time<0,:)=[];
NewData.Time=NewData.Time-NewData.Time(1);

NewData(NewData.Time>1460,:)=[];
writetable(NewData,[subjid '_analyzed.xlsx'],'Sheet','Data');
if(PL_method==2)
    writetable(struct2table(Abs),[subjid '_analyzed.xlsx'],'Sheet','Baseline');
end


NewData2=struct;
NewData2.Time=time_ideal'-Events(1)+300;
NewData2.Event_Rapid=stim(:,2);
NewData2.Event_Hold=stim(:,1);
NewData2.HbO2=HbO2f;
NewData2.Hb=Hbf;
NewData2.CCO=CCOf;
%NewData2.StO2=StO2f;
NewData2.EtCO2=EtCO2;
NewData2.RR=RR;
NewData2.SpO2=SpO2;
NewData2.Pulse=Pulse;
NewData2.Temperature=temp_ideal';
NewData2=struct2table(NewData2);
NewData2(NewData2.Time<0,:)=[];
NewData2.Time=NewData2.Time-NewData2.Time(1);

NewData2(NewData2.Time>1460,:)=[];

writetable(NewData2,[subjid '_analyzed.xlsx'],'Sheet','Data MotionCorrected');
if(PL_method==2)
    writetable(struct2table(Absf),[subjid '_analyzed.xlsx'],'Sheet','Baseline MotionCorrected');
end

fileout=[subjid '_analyzed.pdf'];
import mlreportgen.report.*
import mlreportgen.dom.*
rpt = Report(fileout,'pdf');

tp = TitlePage;
tp.Title = 'FAMPATH Data QC Report';
tp.Author = 'T. Huppert';
add(rpt,tp);
add(rpt,TableOfContents);

ch0 = Chapter;
ch0.Title = 'Info';
sec0=Paragraph;
sec0.WhiteSpace='preserve';
append(sec0,Text(['Subjid: ' subjid ]));
append(sec0,LineBreak);
append(sec0,Text(['Study: ' num2str(info.Study) ]));
append(sec0,LineBreak);
append(sec0,Text(['Date Scanned: ' char(info.Date) ]));
append(sec0,LineBreak);
append(sec0,Text(['Group: ' info.Group ]));
append(sec0,LineBreak);
append(sec0,Text(['Experimenter: ' info.Experimenter ]));
append(sec0,LineBreak);
append(sec0,LineBreak);
append(sec0,LineBreak);
append(sec0,LineBreak);append(sec0,Text(['    NIRS Errors: '  info.NIRS_errors ]));
append(sec0,LineBreak);
append(sec0,Text(['    Eprime Errors'  info.Eprime_errors ]));
append(sec0,LineBreak);
append(sec0,Text(['    Capno Errors: ' info.Capnometry_errors ]));
append(sec0,LineBreak);
append(sec0,Text(['    Other Errors: ' info.Other_errors ]));
append(sec0,LineBreak);
append(sec0,Text(['    Notes: ' info.Behavior_notes ]));
append(sec0,LineBreak);
append(sec0,LineBreak);
append(sec0,Text('Analysis Notes'));

for i=1:length(ProcessingNotes)
    append(sec0,LineBreak);
    append(sec0,Text(ProcessingNotes{i}));
end

append(ch0,sec0);
append(rpt,ch0);





ch1 = Chapter;
ch1.Title = 'Raw Data';

sec1 = Section;
sec1.Title = 'original';

figure;
plot(wavelengths,median(OD_ideal,1),'k','linewidth',2);
set(gca,'Fontsize',14);
xlabel('Wavelength (nm)');
fig(1) = Figure(gcf);
fig(1).Snapshot.Height = '3in';
fig(1).Snapshot.Width = '6in';
fig(1).Snapshot.Caption = 'Raw spectrum (time averaged)';
add(sec1,fig(1));

figure;
plot(time_ideal,median(OD_ideal,2),'k','linewidth',2);
hold on;
ylim=get(gca,'ylim');
ons=find(diff(stim(:,1))>0);
off=find(diff(stim(:,1))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'r','FaceAlpha',.1);
end;
ons=find(diff(stim(:,2))>0);
off=find(diff(stim(:,2))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'b','FaceAlpha',.1);
end;
legend({'Average','Breath','','','','','Rapid'})
set(gca,'Fontsize',14);
set(gca,'Ylim',ylim)
xlabel('Time (s)');
fig(2) = Figure(gcf);
fig(2).Snapshot.Height = '3in';
fig(2).Snapshot.Width = '6in';
fig(2).Snapshot.Caption = 'Raw time course';
add(sec1,fig(2));
add(ch1,sec1);



sec2 = Section;
sec2.Title = 'Corrected';

figure;
plot(wavelengths,median(OD_ideal_filtered,1),'k','linewidth',2);
set(gca,'Fontsize',14);
xlabel('Wavelength (nm)');
fig(3) = Figure(gcf);
fig(3).Snapshot.Height = '3in';
fig(3).Snapshot.Width = '6in';
fig(3).Snapshot.Caption = 'Raw spectrum (time averaged)';
add(sec2,fig(3));

figure;
plot(time_ideal,median(OD_ideal_filtered,2),'k','linewidth',2);
hold on;
ylim=get(gca,'ylim');
ons=find(diff(stim(:,1))>0);
off=find(diff(stim(:,1))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'r','FaceAlpha',.1);
end;
ons=find(diff(stim(:,2))>0);
off=find(diff(stim(:,2))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'b','FaceAlpha',.1);
end;
legend({'Average','Breath','','','','','Rapid'})
set(gca,'Ylim',ylim)
set(gca,'Fontsize',14);
xlabel('Time (s)');
fig(4) = Figure(gcf);
fig(4).Snapshot.Height = '3in';
fig(4).Snapshot.Width = '6in';
fig(4).Snapshot.Caption = 'Raw time course';
add(sec2,fig(4));
add(ch1,sec2);




sec3 = Section;
sec3.Title = 'Physiology';
figure;

subplot(3,1,1);
plot(time_ideal,EtCO2,'k'); xlabel('time (s)'); ylabel('EtCO_2');
hold on;
ylim=get(gca,'ylim');
ons=find(diff(stim(:,1))>0);
off=find(diff(stim(:,1))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'r','FaceAlpha',.1);
end;
ons=find(diff(stim(:,2))>0);
off=find(diff(stim(:,2))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'b','FaceAlpha',.1);
end;
set(gca,'Fontsize',12);
set(gca,'Ylim',ylim)

subplot(3,1,2);
plot(time_ideal,SpO2,'k'); xlabel('time (s)'); ylabel('SpO_2');
hold on;
ylim=get(gca,'ylim');
ons=find(diff(stim(:,1))>0);
off=find(diff(stim(:,1))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'r','FaceAlpha',.1);
end;
ons=find(diff(stim(:,2))>0);
off=find(diff(stim(:,2))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'b','FaceAlpha',.1);
end;
set(gca,'Fontsize',12);
set(gca,'Ylim',ylim)

subplot(3,1,3);
plot(time_ideal,Pulse,'k'); xlabel('time (s)'); ylabel('Pulse Rate');
hold on;
ylim=get(gca,'ylim');
ons=find(diff(stim(:,1))>0);
off=find(diff(stim(:,1))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'r','FaceAlpha',.1);
end;
ons=find(diff(stim(:,2))>0);
off=find(diff(stim(:,2))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'b','FaceAlpha',.1);
end;
set(gca,'Fontsize',12);
set(gca,'Ylim',ylim)
fig(5) = Figure(gcf);
fig(5).Snapshot.Height = '3in';
fig(5).Snapshot.Width = '6in';
fig(5).Snapshot.Caption = 'Raw time course';
add(sec3,fig(5));
add(ch1,sec3);

add(rpt,ch1);

%-----------
ch2 = Chapter;
ch2.Title = 'Chromophore Changes';

sec4 = Section;
sec4.Title = 'Concentration (uncorrected)';

figure;
hold on;
plot(time_ideal,HbO2,'r','linewidth',2)
plot(time_ideal,Hb,'b','linewidth',2)
plot(time_ideal,CCO,'c','linewidth',2)
ylim=get(gca,'ylim');
ons=find(diff(stim(:,1))>0);
off=find(diff(stim(:,1))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'r','FaceAlpha',.1);
end;
ons=find(diff(stim(:,2))>0);
off=find(diff(stim(:,2))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'b','FaceAlpha',.1);
end;

legend({'HbO_2','Hb','CCO','Breath','','','','','Rapid'})

xlabel('time (s)'); ylabel('Concentration [\muM]')
set(gca,'Ylim',ylim)
set(gca,'Fontsize',12);
fig(6) = Figure(gcf);
fig(6).Snapshot.Height = '3in';
fig(6).Snapshot.Width = '6in';
fig(6).Snapshot.Caption = 'Concentration time course';
add(sec4,fig(6));
add(ch2,sec4);


sec5 = Section;
sec5.Title = 'Concentration (corrected)';
figure;
hold on;
plot(time_ideal,HbO2f,'r','linewidth',2)
plot(time_ideal,Hbf,'b','linewidth',2)
plot(time_ideal,CCOf,'c','linewidth',2)
ylim=get(gca,'ylim');
xlabel('time (s)'); ylabel('Concentration [\muM]')
ons=find(diff(stim(:,1))>0);
off=find(diff(stim(:,1))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'r','FaceAlpha',.1);
end;
ons=find(diff(stim(:,2))>0);
off=find(diff(stim(:,2))<0);
for i=1:5;
    patch(time_ideal([ons(i) ons(i) off(i) off(i)]),[ylim(1) ylim(2) ylim(2) ylim(1)],'b','FaceAlpha',.1);
end;
legend({'HbO_2','Hb','CCO','Breath','','','','','Rapid'})
set(gca,'Fontsize',12);
set(gca,'Ylim',ylim)
fig(7) = Figure(gcf);
fig(7).Snapshot.Height = '3in';
fig(7).Snapshot.Width = '6in';
fig(7).Snapshot.Caption = 'Concentration time course';
add(sec5,fig(7));
add(ch2,sec5);

add(rpt,ch2);

if(PL_method==2)
    ch3 = Chapter;
    ch3.Title = 'Baseline';
    sec6 = Section;
    sec6.Title = 'Summary';

    sec6A = Section;
    sec6A.Title = 'Uncorrected';
    tt=summary(struct2table(Abs));
    flds=fields(tt);
    sec6a0=Paragraph;
    sec6a0.WhiteSpace='preserve';

    for i=1:length(flds)
        append(sec6a0,LineBreak);
        append(sec6a0,Text(flds{i}));
        append(sec6a0,LineBreak);
        append(sec6a0,Text(['   Min: ' num2str(tt.(flds{i}).Min)]));
        append(sec6a0,LineBreak);
        append(sec6a0,Text(['   Max: ' num2str(tt.(flds{i}).Max)]));
        append(sec6a0,LineBreak);

        append(sec6a0,Text(['   Median: ' num2str(tt.(flds{i}).Median)]));
        append(sec6a0,LineBreak);
    end
    add(sec6A,sec6a0);

    sec6B = Section;
    sec6B.Title = 'Motion corrected';
    tt=summary(struct2table(Absf));
    flds=fields(tt);
    sec6b0=Paragraph;
    sec6b0.WhiteSpace='preserve';

    for i=1:length(flds)
        append(sec6b0,LineBreak);
        append(sec6b0,Text(flds{i}));
        append(sec6b0,LineBreak);
        append(sec6b0,Text(['   Min: ' num2str(tt.(flds{i}).Min)]));
        append(sec6b0,LineBreak);
        append(sec6b0,Text(['   Max: ' num2str(tt.(flds{i}).Max)]));
        append(sec6b0,LineBreak);

        append(sec6b0,Text(['   Median: ' num2str(tt.(flds{i}).Median)]));
        append(sec6b0,LineBreak);
    end
    add(sec6B,sec6b0);
    add(sec6,sec6A);
    add(sec6,sec6B);
    add(ch3,sec6);
    
    sec7 = Section;
    sec7.Title = 'Time Course';
    sec7A = Section;
    sec7A.Title = 'Uncorrected';
    for i=1:length(flds); Abs.(flds{i})=single(round(Abs.(flds{i})*1000)/1000); Abs.(flds{i})=num2str(Abs.(flds{i}));  end;
    tbl = BaseTable(struct2table(Abs));
    tbl.Title='Uncorrected';
    add(sec7A,tbl);
    add(sec7,sec7A);

    sec7B = Section;
    sec7B.Title = 'Motion corrected';
    for i=1:length(flds); Absf.(flds{i})=single(round(Absf.(flds{i})*1000)/1000); Absf.(flds{i})=num2str(Absf.(flds{i})); end;
    tbl = BaseTable(struct2table(Absf));
    tbl.Title='Motion corrected';
    add(sec7B,tbl);
    add(sec7,sec7B);
    add(ch3,sec7);
    add(rpt,ch3);
end


close(rpt);
rptview(rpt);
close all;

