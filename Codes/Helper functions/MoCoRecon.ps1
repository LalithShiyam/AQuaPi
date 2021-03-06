# This powershell script does the following tasks:
#   1.) Converts the DICOM u-maps with motion to interfile
#   2.) Based on the least temporal difference criteria, the u-maps are assigned to their respective frames.
#   3.) The assigned u-maps are then reconstructed together with their PET counterparts.
#   4.) The decay correction is based on the first time-stamp.
#
# Author - Lalith Kumar Shiyam sundar, M.Sc.
# Quantitative Imaging and Medical Physics, Vienna
# 20.02.2018
#-------------------------------------------------------------------------------------------------------------------#
#                                                       PROGRAM START
#-------------------------------------------------------------------------------------------------------------------#

# Initialization zone
#---------------------#
function MoCoRecon
{
Param (
$subjectID, # this is volunteer code
$processingFolder, # Usually this is going to be the recons desktop, where the pet raw data and the motion reflected attenuation maps are stored.
$attenuationMapToCopy,# name of the attenution map (which contains the patient motion).
$pathToDynUmap, # path of the dynamic attenuation maps.
$pathOfListModePET, # Path of the list-mode raw data.
$pathOfParamFile, # path of the parameter file for reconstruction.
$pathOfPETMRlink, # path of the PET-MR navigator correspondence text file.
$pathToCopyReconPETdata, # path to copy the recon data.
[int]$samePatient
)

# Hard-coded variables  
#---------------------# 
$WhereAmI='C:\Users\admin\Desktop\' # Keep it hardcoded
$motherFolder=Join-Path -path $WhereAmI -ChildPath $subjectID # This is the path where the processing will occur
$workingDir='e7_WorkDir' # folder which contains the folder PET list-mode data and "AC"
$AttenuationDir='AC' # folder which will be created to store one attenuation map at a time, just done to trick the siemen's reconstruction software.
$motionUmapsDirIF='Interfile_Motion_uMaps' # folder which contains the interfile version of the dynamic attenuation maps.
$motionUmapsDirDicom='Dicom_Motion_uMaps' # folder which contains the dynamic ac maps based on the mr navigator's motion vectors.
$reconstructedPETData='MoCo_Recons_PET_' + $attenuationMapToCopy # folder which contains the reconstructed PET data.
$samePatientTrigger=0;

# Hard-coding stops here, paths created based on the previous inputs.
#-------------------------------------------------------------------#

$workingPath=Join-Path -path $motherFolder -ChildPath $workingDir
$attenuationPath=Join-Path -path $workingPath -ChildPath $AttenuationDir # path containing the AC map in the working directory
$pathOfMotionUmaps= Join-Path -path $motherFolder -ChildPath $motionUmapsDirIF # path of the motion u-maps folder in interfile format.
$pathOfReconsPET=Join-Path -Path $motherFolder -ChildPath $reconstructedPETData # path of the reconstructed PET inside the working directory of the subject.
$pathOfDicomMotionuMaps=Join-Path -path $motherFolder -ChildPath $motionUmapsDirDicom   # path of the motion u-maps folder in dicom format.
if ($samePatient -eq $samePatientTrigger)
    {
    New-Item -ItemType directory -Path $motherFolder
    }
New-Item -ItemType directory -path $pathOfDicomMotionuMaps    
New-Item -ItemType directory -path $pathOfMotionUmaps # directory where the motion filled interfile u-maps are stored
if ($samePatient -eq $samePatientTrigger)
    {
    New-Item -ItemType directory -path $workingPath # creates e7_WorkDir folder inside mother folder 
    }
New-Item -ItemType directory -path $attenuationpath # creates AC folder inside the e7_WorkDir
New-Item -ItemType directory -Path $pathOfReconsPET # creates folder for reconstructed pet data
Copy-Item -Force -Recurse -Verbose $pathToDynUmap\* -Destination $pathOfDicomMotionuMaps   # copies the dynamic attenuation maps to the DICOM motion maps folder in the mother folder.
if ($samePatient -eq $samePatientTrigger)
   {
    Copy-Item -Force -Recurse -Verbose $pathOfListModePET\* -Destination $workingPath # copies the list-mode data to the e7_WorkDir
   }
Copy-Item -Force -Recurse -Verbose $pathOfParamFile -Destination $workingPath # copies the parameter file to the e7_WorkDir
Copy-Item -Force -Recurse -Verbose $pathOfPETMRlink -Destination $workingPath # copies the pet-mr navigator correspondence text file to the mother folder



# Count the number of Dicom_Motion_uMaps available and start the JSrecon to convert the DICOM u-maps to interfile u-maps
#-----------------------------------------------------------------------------------------------------------------------#

cd($WhereAmI)
$count=0

foreach ($i in Get-ChildItem -path $pathOfDicomMotionuMaps)
{
    $numOfFolders=(dir $pathOfDicomMotionuMaps | measure).Count;
    $uMapOfInterest=Join-Path -Path $pathOfDicomMotionuMaps -ChildPath $i.Name
    Write-Host "Copying the right u-map to the working directory..."
    Copy-Item -Force -Recurse -Verbose $uMapOfInterest\* -destination $attenuationPath # copies the dicom u-maps to the e7_WorkDir
    Write-Host "Running JSrecon12..."
    cscript.exe C:\JSRecon12\JSRecon12.js $workingPath $pathOfParamFile  #passing the parameters to JSrecon12.
    $uMapseriesName= -Join("IF_",$i) # creates a folder name for the converted u-maps in interfile format, index starts with '00'
    $IFmotionSeriesPath= Join-Path -path $pathOfMotionUmaps -ChildPath $uMapseriesName # creates a path to store the interfile
    #New-Item -ItemType directory -Path $IFmotionSeriesPath
    $convertedFolder=Get-ChildItem -path $motherFolder -filter *converted* # finding the e7_WorkDir_Converted folder
    $pathOfConvertedFolder=Join-Path -path $motherFolder -ChildPath $convertedFolder # creating a path for the converted folder 
    $interfileUMapFolder=Get-ChildItem -path $pathOfConvertedFolder -filter *LM-00* # Looking for the directory which contains the interfile u-maps
    $pathToIFuMapFolder=Join-Path -path $pathOfConvertedFolder -ChildPath $interfileUMapFolder # preparing the folder path for initiating transfer
    Copy-Item -Force -Recurse -Verbose $pathToIFuMapFolder -Filter *umap* -Destination $IFmotionSeriesPath # transfers the interfile u-maps from working directory to the motion series directory
    $oldACfileName=Get-ChildItem -Path $IFmotionSeriesPath -Filter *.v # find the old attenuation map names # getting the old u-map names which need to be replaced by new names in the main and local header.
    Get-ChildItem -Path $IFmotionSeriesPath -Recurse |Rename-Item -NewName{$_.name -replace $workingDir,$uMapseriesName} # renaming the old u-map names with new u-map names.
    $newACfileName=Get-ChildItem -Path $IFmotionSeriesPath -Filter *.v # find the new attenuation map names 
    $mainHdrName=Get-ChildItem -Path $IFmotionSeriesPath -Filter *.mhdr #find the main header. # includes the hardware attenuation maps and human u-maps
    $localHdrName=Get-ChildItem -Path $IFmotionSeriesPath -filter *.hdr # find the local header # same as above
    $temp=0
    Write-Host "Rewriting the file names inside the old header with new ones..."
    foreach ($j in $mainHdrName) # changing the old u-map names in main header and local header of the hardware/human u-map.
    {
    $tempPath=Join-Path -path $IFmotionSeriesPath -ChildPath '\'
    $mainHdrPath=Join-Path -Path $tempPath -ChildPath $mainHdrName[$temp].Name
    (Get-Content $mainHdrPath) -replace ($oldACfileName[$temp].Name,$newACfileName[$temp].Name) | Set-Content $mainHdrPath
    $localHdrPath=Join-Path -Path $tempPath -ChildPath $localHdrName[$temp].Name
    (Get-Content $localHdrPath) -replace ($oldACfileName[$temp].Name,$newACfileName[$temp].Name) | Set-Content $localHdrPath 
    $temp=$temp+1;
   } 
   $uMapfiles2Remove=Join-Path -path $attenuationPath -ChildPath '\*'
   Remove-Item -Path  $uMapfiles2Remove | Where {! $_.PSIsContainer}  # old u-map files are removed from the AC folder of working directory
   Write-Host "Removing the old u-map files..."
   $count=$count+1
   if ($count -eq $numOfFolders ) 
    {
      Write-Host "Not removing the converted folder..."
      Write-Host $numOfFolder
    }
      else
      {
        Remove-Item -Path $pathOfConvertedFolder -recurse 
        Write-Host "Removing the converted folder..." 
      } 
   
 }   
    
# Perform the histogramming of the "all-sinogram" by running the histogram batch file
#-------------------------------------------------------------------------------------------
$histogrammingBatchFile=Get-ChildItem -Path $pathToIFuMapFolder -Filter *Histogram* # searching for the bat file which contains the keyword "histogram"
$pathToHistBatFile=Join-Path -Path $pathToIFuMapFolder -childPath $histogrammingBatchFile.Name # creating the path for the histogramming bat file.
cd $pathToIFuMapFolder # Must be in the path to store the output in this particular folder
Start-Process $pathToHistBatFile -Wait  # running the histogramming bat file.
Write-Host Data stored in $pathToIFuMapFolder !

# Code snippet for opening the text file and modifying the PET sinogram header.
#--------------------------------------------------------------------------------

cd $pathToIFuMapFolder
$sinoMainHdrName=Get-ChildItem -Path $pathToIFuMapFolder -Filter *sino.mhdr 
$sinoMhdrPath=Join-Path $pathToIFuMapFolder -ChildPath $sinoMainHdrName

# Read in the PET-MR correspondence file
#-----------------------------------------

$MRlink=Get-Content($pathOfPETMRlink)

# Create a function for changing the sinogram headers
#--------------------------------------------------------------------------------

function Get-changeMainHeaders 
{
    Param (
    $numOfPETframes,
    $startPETframe, 
    $endPETframe, 
    $sinoMhdrPath,
    $sinogramFiles
    )
 
 # indexing of powershell starts from zero, tiny modification. 
 

 $startPETframe=$startPETframe
 $endPETframe=$endPETframe
 
 # part where the total number of datasets are changed
 
 $stringToReplace_1='number of time frames:=.+'
 $stringToReplace_2='!total number of data sets:=.+'
 
 $modifiedString_1='number of time frames:=' + $numOfPETframes
 $modifiedString_2='!total number of data sets:=' + $numOfPETframes
 
(Get-Content ($sinoMhdrPath)) | ForEach-Object{$_ -replace $stringToReplace_1, 
$modifiedString_1} | Set-Content ($sinoMhdrPath)
(Get-Content ($sinoMhdrPath)) | ForEach-Object{$_ -replace $stringToReplace_2, 
$modifiedString_2} | Set-Content ($sinoMhdrPath)

# part where the 'dataset[n]' region will be added according to the needs of the 
# program 

<#$searchString=[regex]::Escape("%data set [")
$sinogramFiles=Get-Content $sinoMhdrPath |Select-String -Pattern $searchString | Select-Object -ExpandProperty 'Line' 
$sinogramFilesIndex=Get-Content $sinoMhdrPath |Select-String -Pattern $searchString | Select-Object -ExpandProperty 'LineNumber' #>

# Removing and adding the sinogram header lines present at the end of the sino.mhdr 

Set-Content -Path $sinoMhdrPath -Value (get-content -Path $sinoMhdrPath | Select-String -Pattern $searchString -NotMatch) #remove 
#the multiple datasets at the end of the text file, to replace it with my own files.

for ($i=1; $i -le $numOfPETframes; $i++)
{
$strCount=[string]$i
$stringToAdd="%data set ["+$strCount+"]:"
$tempFile=$sinogramFiles[$startPETframe] -replace ".+:",$stringToAdd
Add-Content -path $sinoMhdrPath $tempFile
$startPETframe=$startPETframe+1 
}
}


# Loop over and find out if the row value changes.
#-------------------------------------------------------

# Initialization for iteration.

$startIndex=0
$nextIndex=1
$TotalNumOfPETframes= Get-Content $pathOfPETMRlink | Measure-Object -Line
$startPETframe=0
$numOfPETframes=1
$searchString=[regex]::Escape("%data set [")
$sinogramFiles=Get-Content $sinoMhdrPath |Select-String -Pattern $searchString | Select-Object -ExpandProperty 'Line' 
$sinogramFilesIndex=Get-Content $sinoMhdrPath |Select-String -Pattern $searchString | Select-Object -ExpandProperty 'LineNumber' 

for ($i=1; $i -le $TotalNumOfPETframes.Lines; $i++)
{
    if ($MRlink[$startIndex] -eq $MRlink[$nextIndex])
    {
        $startIndex++
        $nextIndex++
    }
    else
    {
        $endPETframe=$startIndex
        $numOfPETframes= ($endPETframe-$startPETframe)+1
        $associatedMRnav=$MRlink[$startIndex]
        Get-changeMainheaders $numOfPETframes $startPETframe $endPETframe $sinoMhdrPath $sinogramFiles # <- modifying the sinogram headers 
        $oldUmapfiles2Remove=Join-Path -path $pathToIFuMapFolder -ChildPath '\*umap*'  
        Remove-Item -Path  $oldUmapfiles2Remove | Where {! $_.PSIsContainer}  # old CT files are removed from the AC folder of working directory
        Write-Host "Removing the old u-maps..."
        $uMapStringToSearch='*_'+ $associatedMRnav
        $UMapFolderOfInterest=Get-ChildItem -path $pathOfMotionUmaps -filter $uMapStringToSearch # Looking for the directory which contains the interfile u-maps
        $pathToUMapFolderOfInterest=Join-Path -path $pathOfMotionUmaps -ChildPath $UMapFolderOfInterest\* 
        Copy-Item -Path $pathToUMapFolderOfInterest -Destination $pathToIFuMapFolder
        $uMapMainHdrName=Get-ChildItem -Path $pathToUMapFolderOfInterest -Filter *umap.mhdr
        $uMapHardwareHdrName=Get-ChildItem -Path $pathToUMapFolderOfInterest -Filter *hardware.mhdr
        $stringToReplace='-u.+'
        $modifiedString= '-u "'+$uMapMainHdrName.Name+'"'+','+'"'+$uMapHardwareHdrName.Name+'"'
        $reconFileName=Get-ChildItem -Path $pathToIFuMapFolder -Filter *OP.bat 
        $pathToReconBat=Join-Path $pathToIFuMapFolder -ChildPath $reconFileName
        (Get-Content ($pathToReconBat)) | ForEach-Object{$_ -replace $stringToReplace, 
        $modifiedString} | Set-Content ($pathToReconBat)
        # Perform the reconstruction
        #-------------------------------------------------------------------------------------------
        cd $pathToIFuMapFolder # Must be in the path to store the output in this particular folder
        Start-Process -FilePath $pathToReconBat -Wait # running the OP bat file.
        $IF2DcmBat=Get-ChildItem -Path $pathToIFuMapFolder -Filter *IF2Dicom.bat
        $pathToIF2DcmBat=Join-Path $pathToIFuMapFolder -ChildPath $IF2DcmBat
        Start-Process -FilePath $pathToIF2DcmBat -Wait 
        Write-Host Data stored in $pathToIFuMapFolder !
        Write-Host $associatedMRnav $numOfPETframes $startPETframe $endPETframe
        $oldDicomFolderPath=Get-ChildItem -Recurse | ?{ $_.PSIsContainer } | % { $_.FullName }
        $reconstructedPET='PET-Frame-'+$startPETframe+'-'+$endPETframe
        Rename-Item $oldDicomFolderPath $reconstructedPET
        Move-Item $reconstructedPET $pathOfReconsPET
        $startIndex=$endPETframe+1
        $nextIndex=$startIndex+1
        $startPETframe=$startIndex
     }
}
    

copy-Item -Force -Recurse -Verbose $pathOfReconsPET $pathToCopyReconPETdata
Get-ChildItem $pathOfReconsPET -recurse -force -verbose .| Remove-Item -force -recurse -verbose
Remove-Item -Force -Recurse -Verbose $pathOfReconsPET
Get-ChildItem $attenuationPath -recurse -force -verbose .| Remove-Item -force -recurse -verbose
Remove-Item -Force -Recurse -Verbose $attenuationPath 
Get-ChildItem $pathOfMotionUmaps -recurse -force -verbose .| Remove-Item -force -recurse -verbose
Remove-Item -Force -Recurse -Verbose $pathOfMotionUmaps
Get-ChildItem $pathOfDicomMotionuMaps -recurse -force -verbose .| Remove-Item -force -recurse -verbose
Remove-Item -Force -Recurse -Verbose $pathOfDicomMotionuMaps
Get-ChildItem $pathOfConvertedFolder -recurse -force -verbose .| Remove-Item -force -recurse -verbose
Remove-Item -Force -Recurse -Verbose $pathOfConvertedFolder
Remove-Item -Force -Recurse -Verbose $workingPath\*.txt*
}










    