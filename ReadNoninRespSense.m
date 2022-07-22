
EtCO2=[];
CO2=[];
RR=[];

startIndices=find(data==hex2dec('FF'));
for i=1:length(startIndices)
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
    
%     checksum=256*local_data(5)+local_data(6);
%     checksum2=local_data(3)+local_data(4);
%     if(checksum~=checksum2)
%         error('checksum missmatch');
%     end
    
end