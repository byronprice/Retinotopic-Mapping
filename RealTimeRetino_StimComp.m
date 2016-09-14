% RealTimeRetino_StimComp

startEXP = 254;
endEXP = 255;

startRUN = 252;
endRUN = 253;

numChans = 2;

tcpipClient = tcpip('128.197.59.169',30000,'NetworkRole','client');
bufferSize = 50000; % bytes, a big number (won't need this much)
set(tcpipClient,'InputBufferSize',bufferSize);
set(tcpipClient,'Timeout',5);
fopen(tcpipClient);

WaitSecs(20);
usb = usb1208FSPlusClass;
usb.strobeEventWord(startEXP);
bigData = [];
for ii=1:numChans
    for jj=1:10
        usb.strobeEventWord(jj);
        dataSize = fread(tcpipClient,3,'double');
        
        if isempty(dataSize) == 0
            temp = fread(tcpipClient,dataSize(1)*dataSize(2),'double');
            data = reshape(temp,[dataSize(1),dataSize(2)]);
            bigData = [bigData;data];
        end
    end
    usb.strobeEventWord(endRUN);
end

bigData