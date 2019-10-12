function [PETframingMidPoint]=getPETmidTime()

PETframing=[5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,60,120,300,300,300,300,300,300,300,300,300,300];
PETframing=[0,PETframing,300];

% creating the mid-time pet framing. 

time=cumsum(PETframing); time=time';
timeShiftOne=circshift(time,-1);
timeAvg=(time+timeShiftOne)./2;
timeAvg=timeAvg(1:end-1);
PETframingMidPoint=timeAvg;
end