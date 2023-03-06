load FAMPATH_REF_CAL.mat

Darkf=Dark-ones(size(Dark,1),1)*median(Dark(1:200,:),1);
Darkf=medfilt2(Darkf,[5 5]);
Dark_noise=median(Darkf,2);

RoomOnf=RoomOn-ones(size(RoomOn,1),1)*median(RoomOn(1:200,:),1);
RoomOnf=medfilt2(RoomOnf,[5 5]);
RoomOn_noise=median(RoomOnf,2);

RoomOfff=RoomOff-ones(size(RoomOff,1),1)*median(RoomOff(1:200,:),1);
RoomOfff=medfilt2(RoomOfff,[5 5]);
RoomOff_noise=median(RoomOfff,2);


for i=1:85; 
    tmp=LightOff_calib(:,:,i)-ones(size(LightOff_calib,1),1)*median(LightOff_calib(1:200,:,i),1);
    tmp=medfilt2(tmp,[5 5]);
    A(:,i)=median(tmp,2)-DarkNoise; 
end;