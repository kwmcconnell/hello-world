
% Automated R-Mode Processing of Native Data Echocardiograms
%{
Kelly McConnell December 2017

The following program reads native data DICOM echocardiogram files to
reorder selective images in respiration mode. The program locates end
systole markers and R wave peak markers. The user can input a shift, which
the program uses to find images that distance away from the R wave peak. The
program then reorders the images based on the respiration waveform, inhale to
exhale. The program produces 3 outputs: a plot of both the ECG and respiration
waveforms, a video of the full 2D echo, and a video of the R Mode echo
%}
%{
NOTE: How to select correct DICOM file(s) for processing.
1) Open desired study in Xcelera.
2) Note which loop(s) you want to process (it needs to have "raw data"
yellow triangle and both ECG and respiratory waveforms).
3) In the Xclera Study Viewer, right click on study with correct date and
select "Copy Study" - DICOM will be copied to C:\LocalStorage directory
4) Open the newly created folder within LocalStorage directory
5) Look for files with large sizes (due to including raw data). 
6) Double-click on the filename (selecting DICOM viewer/Showcase as
program) to make sure this is the desired loop.
7) Copy this file into MATLAB folder (note the MATLAB program will ask for the
filename)
%}

close all %close all open figures
clear all %clear workspace
%Input parameters
prompt = 'Input image name and press enter (Example: 6OSUDFF3): ';
imagename = input(prompt, 's'); %'6OSUDFF3'; %DICOM image file
err = 10000; %+- between data timestamp and image timestamp (image period is about 20,000)
SystoleShift = 0; %shift in number of frames from End Systole marker
RWaveShift = .4; %shift in percent of R wave period from R wave marker
percentSmooth = .1; %how much of the wave the smoothing algorithm works over

%program reads in data from DICOM
try
    fileinfo = dicominfo(imagename); %file information
catch
    warning('Image name not found. make sure the image is in the Matlab folder and it is spelled correctly. Please rerun the file');
end
[data,params, reserved, datapos] = pitDicomNativeExtract(imagename, 'UDM_USD_DATATYPE_DIN_PHYSIO'); %read in physio data
[data2,params2, reserved2, datapos2] = pitDicomNativeExtract(imagename, 'UDM_USD_DATATYPE_DIN_2D_ECHO'); %read in 2D echo image data
imagedata = dicomread(imagename); %read in 2D echo images

%set up echo image data
sizeImageData = size(imagedata);
numImages = sizeImageData(4); %number of image frames
time_images =datapos2;% %timestamps for echo images
dataPointsArray = ECGDataProcessing(imagename,data, params, RWaveShift, SystoleShift, true);

%chosing shift
prompt = 'Pick the shift to use: r wave or systole (case sensitive): ';
shiftChoice = input(prompt, 's');
if isempty(strfind(shiftChoice, 'r wave')) == false
    useRWaveShift = true;
elseif isempty(strfind(shiftChoice, 'systole')) == false
        useRWaveShift = false;
end
done = false;
while done == false
    prompt = ['input shift value (current values: R wave = ' num2str(RWaveShift) ', Systole = ' num2str(SystoleShift) ')' \n ];
    shiftValue = input(prompt);
    if useRWaveShift == true
        RWaveShift = shiftValue;
    else
        SystoleShift = shiftValue;
    end
    dataPointsArray = ECGDataProcessing(imagename,data, params, RWaveShift,SystoleShift, useRWaveShift);
    prompt = 'Do you want to use this shift value? please enter 1 for yes or 0 for no';
    isDone = input(prompt);
    if isDone == '1'
        done = true;
    end
end
done = false;
RespDataProcessing(dataPointsArray, percentSmooth)
while done == false
    prompt = ['input smooth value (current value:  ' num2str(percentSmooth) ')' \n ];
    percentSmooth = input(prompt);
    RespDataProcessing(dataPointsArray, percentSmooth, err, time_images)
    prompt = 'Do you want to use this smooth value? please enter 1 for yes or 0 for no';
    isDone = input(prompt);
    if isDone == '1'
        done = true;
    end
end
    


function [dataPointsArray] = ECGDataProcessing(imagename,data, params, RWaveShift,SystoleShift, useRWaveShift)

%set up physio data
numFrames = double(params(2)); %number of frames of ecg/respiratory data
xxyy = squeeze(data); %remove singleton dimentions
time_data = (xxyy(1,:));%/(1000*fileinfo.CineRate); %timestamps for physio data
ecg = xxyy(5,:); %ecg wave
resp = xxyy(6,:); %respiratory wave
if size(resp) == size(find(resp == 0))
    error('respiratory data not found');
end


DataFormatArray = dec2bin(data(3,:)); %data format information converted to binary
EndSystoleArray = zeros(1,numFrames);
RWaveArray = zeros(1,numFrames);


%index arrays
EndSystoleArray_shift = [0,0,0,0,0,0,0]; %[index, time, ecg value, resp value, resp', resp'', resp smooth]
RWaveArray_shift = [0,0,0,0,0,0,0]; %[index, time, ecg value, resp value, resp', resp'', resp smooth]
resp_index = zeros(1,numFrames);

%finding End Systole points and R wave peaks
n = 1; %systole index
m = 1; %r wave index
for i = 1:numFrames
    
    if  DataFormatArray(i,:) == '10000100000000' %indicates End Systole point
        
        EndSystoleArray(n,:) = [i,time_data(i), ecg(i), resp(i), resp_smooth_dt(i),resp_smooth_dtdt(i),resp_smooth(i)];
        n = n + 1;
    elseif DataFormatArray(i,:) == '00000101000000' %indicates r wave peak point
        RWaveArray(m,:) = [i,time_data(i), ecg(i), resp(i), resp_smooth_dt(i),resp_smooth_dtdt(i), resp_smooth(i)];
        m = m+1;
    end
end
%getting images at SystoleShift% of cardiac cycle
index1 = EndSystoleArray(1,1);
index2 = EndSystoleArray(2,1); %absolute shift
for i = 1:n-1
    s = round(EndSystoleArray(i,1) + (index2-index1)*SystoleShift);
    EndSystoleArray_shift(n,:) = [s,time_data(s), ecg(s), resp(s), resp_smooth_dt(s),resp_smooth_dtdt(s),resp_smooth(s)];
end

%getting images at RWaveShift% of R Wave
index1 = RWaveArray(1,1);
index2 = RWaveArray(2,1); %absolute shift
for i = 1:m-1
    z = round(RWaveArray(i,1) + (index2-index1)*RWaveShift);
    RWaveArray_shift(i,:) = [z, time_data(z), ecg(z), resp(z), resp_smooth_dt(z),resp_smooth_dtdt(z), resp_smooth(z)];
end

if useRWaveShift == true 
    dataPointsArray = RWaveArray_shift;
else
    dataPointsArray = EndSystoleArray_shift;
end



%% Plots

%ECG Wave plots
figure()
figure('Units','inches', 'Position', [1 1 10 6], 'PaperPositionMode', 'auto');
plot(time_data,ecg) %ecg plot
axis tight
ylabel('ECG wave')
title([imagename, ' Native Data'])
hold on
%add end systole and r wave peak points to ecg plot
plot(time_data(EndSystoleArray(:,1)), EndSystoleArray(EndSystoleArray(:,1)), 'o') %end systole points
hold on
% plot(time_data(EndSystoleArray_index(:,1)+SystoleShift), EndSystoleArray(EndSystoleArray_index(:,1)+SystoleShift), 'o') %end systole points + shift
% hold on
plot(time_data(RWaveArray(:,1)), RWaveArray(RWaveArray(:,1)), '*') %r wave points
hold on
plot(time_data(dataPointsArray(:,1)), ecg(dataPointsArray(:,1)), '.', 'MarkerSize', 20) %image sample points
hold on
legend('ecg data', 'End Systole Points', 'R Wave Peaks', 'Location','eastoutside', ['Image Sample (R wave percent shift: ' num2str(RWaveShift) ', End systole shift: ' num2str(SystoleShift)])


end
%%%need to figure out smoothing stuff, end 4 points in data point arrays
%%%are from the smoothing, either new array or add to it here
function [] = RespDataProcessing(dataPointsArray, percentSmooth, err, time_images)
xxyy = squeeze(data); %remove singleton dimentions
time_data = (xxyy(1,:));%/(1000*fileinfo.CineRate); %timestamps for physio data
resp = xxyy(6,:); %respiratory wave
%smoothing resp wave
resp_smooth = smooth(time_data, resp, percentSmooth, 'lowess'); %lowess method- regression smoothing over 10% of wave
resp_smooth_dt = diff(resp_smooth); %finding slope of the respiratory wave
resp_smooth_dtdt = diff(resp_smooth,2); %finding 2nd derivative of respiratory wave

%finding max/min for resp data
[peaks,loc] = findpeaks(resp_smooth); %finding maximums of respiratory wave
[peaks_n, loc_n] = findpeaks(-resp_smooth); %finding minimums of respiratory wave
[peaks_dt,loc_dt] = findpeaks(resp_smooth_dt); %finding maximuns of respiratory slope

%labels for resp plot (A, B...)
labels = {};
for i = 1:m-1 %%n for systole, m for r wave
    labels{i} = char(64+i); %labels starting at A
end

%determining where the points are in the resp cycle
percentResp = zeros(1,m-1);
loc_points = [loc;loc_n]; %put all the max and min points in one list
[loc_points_sort, i_loc_points_sort]= sort(loc_points); % put them in order
if n<m %make sure index does not go out of bounds if using systole array
    m = n;
end

for i = 1:m-1
    
    ind = dataPointsArray(i,1); 
    done = false;
    j = 1;
    while done == false %finding what two points the image is between (point1<ind<point2)
        if j > max(size(loc_points_sort)) %end boarder case
            point1 = loc_points_sort(j-1);
            point2 = numFrames;
            done = true;
        else
            if ind < loc_points_sort(j)
                done = true;
                point2 = loc_points_sort(j);
                if j == 1 %beginning boarder case
                    point1 = 1;
                else
                    point1 = loc_points_sort(j-1);
                end
                
            end
        end
        j = j+1;
    end
    percentResp(i) = (ind-point1)/(point2-point1);  
end
j = 1; %inhale index
k = 1; %exhale index

for i = 1:m-1
    slope = dataPointsArray(i,5);
    if slope > 0
        inhale(j) = [percentResp(i)]; %finding points on inhale
        j = j+1;
    else
        exhale(k) = [percentResp(i)]; %finding points on exhale
        k = k+1;
    end
end
%ordering the respiratory cycle based on percentResp
[inhale_sort,i_inhale] = sort(inhale);
[exhale_sort,i_exhale] = sort(exhale);
resp_order = [inhale_sort,exhale_sort]; %ordered resp points from inhale-->exhale

%ordered index and label lists
for i = 1:m-1
    resp_order_index(i) = find(resp_order(i)==percentResp);
    resp_order_label(i) = labels{find(resp_order(i)==percentResp)};
end
resp_order_label
%creating ordered timestamp array
for i = 1:m-1
    resp_order_timestamps(i) = dataPointsArray(resp_order_index(i),2);  %%%change between RWaveArray_index and EndSystoleArray_index
end
% echo image ordering
for i = 1:m-1
    time_marker= dataPointsArray(resp_order_index(i),2); %%%change between RWaveArray_index and EndSystoleArray_index
    [i_echo] = find(abs(time_marker-time_images)<err, 1, 'first');
    image_time(i) = time_images(i_echo);
    image_index(i) = i_echo;
end
figure()
figure('Units','inches', 'Position', [1 1 10 6], 'PaperPositionMode', 'auto');
plot(time_data,resp) %plot respiratory wave
axis tight
ylabel('Respiratory Wave')
xlabel('Timestamp')
title([imagename, ' Native Data'])
hold on
plot(time_data, resp_smooth) %plot smoothed respiratory wave
hold on
plot(time_data(loc), resp_smooth(loc), '.', 'MarkerSize', 20) %plot respiratory maximum
hold on
plot(time_data(loc_n), resp_smooth(loc_n), '.', 'MarkerSize', 20) %plot respiratory minimum
hold on
%add end systole points to resp graph
plot(time_data(dataPointsArray(:,1)), dataPointsArray(:,7), '*')
hold on
text(time_data(dataPointsArray(:,1)), dataPointsArray(:,7), labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom') %label end systole points
legend('respiratory data','smoothed respiratory data', 'maximums', 'minimums', 'image marker', 'Location','eastoutside')
end

