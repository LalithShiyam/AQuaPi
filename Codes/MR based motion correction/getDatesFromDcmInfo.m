function [dateVector]=getDatesFromDcmInfo(pathOfReconPET)
cd(pathOfReconPET)
fNames=dir('*.ima')
if isempty(fNames)
    fNames=dir('*.dcm')
end
DCMinfo=dicominfo(fNames(1).name);
studyDate(lp,:)=[DCMinfo.StudyDate(7:8),'.',DCMinfo.StudyDate(5:6),'.',DCMinfo.StudyDate(1:4)];
dateVector(lp,:)=datetime(studyDate(lp,:),'InputFormat','dd.MM.yyyy');
end