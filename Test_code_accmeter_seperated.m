clear;

% Adds directories to MATLAB path

% Paths for the linear method
addpath(fullfile(pwd, 'Linear'));
addpath(fullfile(pwd, 'Util'));
addpath(fullfile(pwd, 'matlabPyrTools'));
addpath(fullfile(pwd, 'matlabPyrTools', 'MEX'));

% Paths for the phase-based method
addpath(fullfile(pwd, 'PhaseBased'));
addpath(fullfile(pwd, 'pyrToolsExt'));
addpath(fullfile(pwd, 'Filters'));

%%
dataDir = './Detected_motions';

resultsDir = 'ResultsSIGGRAPH2013/';
mkdir(resultsDir);
% defaultPyrType = 'octave';
% defaultPyrType = 'halfOctave'; % Half octave pyramid is default as discussed in paper
% defaultPyrType = 'smoothHalfOctave';
defaultPyrType = 'quarterOctave';
%% Crane
inFile = fullfile(dataDir, 'VID_20191215_082949.mp4_1st_122_127.avi');
samplingRate = 60; % Hz
loCutoff = 2;    % Hz
hiCutoff = 3;    % Hz
alpha = 10;    
sigma = 3;         % Pixels
temporalFilter = @FIRWindowBP; 
pyrType = defaultPyrType;

vidFile=inFile;
magPhase=alpha;
fl=loCutoff;
fh=hiCutoff;
fs=samplingRate;
outDir=resultsDir;
varargin = {'sigma', sigma,'pyrType', pyrType,'temporalFilter', temporalFilter,'scaleVideo', 1};

%% Read Video
    vr = VideoReader(vidFile);
%% crop video
                    vidObj = VideoWriter('testvideo.avi');
                    vidObj.Quality = 100;
                    vidObj.FrameRate = vr.FrameRate;
                    open(vidObj);
                    c = read(vr, 1);
                    [~,RECT] = imcrop(c);% cropping the frame [XMIN YMIN WIDTH HEIGHT]
                    close(figure(1));
                    for img =1:vr.NumberOfFrames
                        b = read(vr, img);
                        b2 = imcrop(b,RECT);% cropping the frame [XMIN YMIN WIDTH HEIGHT]
                        writeVideo(vidObj,b2)

                    end
                    close(vidObj);
                    vr = VideoReader('testvideo.avi');
%%
                    
    [~, writeTag, ~] = fileparts(vidFile); % return main name of the file('crane')
    FrameRate = vr.FrameRate;    
    vid = vr.read();
    [h, w, nC, nF] = size(vid); %nC=Number of colors, nF=Number of frames
    
%% Parse Input
    p = inputParser();

    defaultAttenuateOtherFrequencies = false; %If true, use reference frame phases
    pyrTypes = {'octave', 'halfOctave', 'smoothHalfOctave', 'quarterOctave'}; 
    checkPyrType = @(x) find(ismember(x, pyrTypes));
    defaultPyrType = 'octave';
    defaultSigma = 0;
    defaultTemporalFilter = @FIRWindowBP;
    defaultScale = 1;
    defaultFrames = [1, nF];
    
    addOptional(p, 'attenuateOtherFreq', defaultAttenuateOtherFrequencies, @islogical);
    addOptional(p, 'pyrType', defaultPyrType, checkPyrType);
    addOptional(p, 'sigma', defaultSigma, @isnumeric);   
    addOptional(p, 'temporalFilter', defaultTemporalFilter);
    addOptional(p, 'scaleVideo', defaultScale);
    addOptional(p, 'useFrames', defaultFrames);
    
    parse(p, varargin{:});

    refFrame = 1;
    attenuateOtherFreq = p.Results.attenuateOtherFreq;
    pyrType            = p.Results.pyrType;
    sigma              = p.Results.sigma;
    temporalFilter     = p.Results.temporalFilter;
    scaleVideo         = p.Results.scaleVideo;
    frames             = p.Results.useFrames;

     %% Compute spatial filters        
    vid = vid(:,:,:,frames(1):frames(2)); %from 1st frame to last frame (frames(2)=last frame)
    [h, w, nC, nF] = size(vid);
    if (scaleVideo~= 1) %Scales the video
        [h,w] = size(imresize(vid(:,:,1,1), scaleVideo));
    end
    
    
    fprintf('Computing spatial filters\n');
    ht = maxSCFpyrHt(zeros(h,w)); % maxSCFpyrHt=Specifies the maximum number of octaves that can be in a steerable pyramid of image IM.
    switch pyrType
        case 'octave'
            filters = getFilters([h w], 2.^[0:-1:-ht], 4);
            repString = 'octave';
            fprintf('Using octave bandwidth pyramid\n');        
        case 'halfOctave'            
            filters = getFilters([h w], 2.^[0:-0.5:-ht], 8,'twidth', 0.75);
            repString = 'halfOctave';
            fprintf('Using half octave bandwidth pyramid\n'); 
        case 'smoothHalfOctave'
            filters = getFiltersSmoothWindow([h w], 8, 'filtersPerOctave', 2);           
            repString = 'smoothHalfOctave';
            fprintf('Using half octave pyramid with smooth window.\n');
        case 'quarterOctave'
            filters = getFiltersSmoothWindow([h w], 8, 'filtersPerOctave', 4);
            repString = 'quarterOctave';
            fprintf('Using quarter octave pyramid.\n');
        otherwise 
            error('Invalid Filter Types');
    end

    [croppedFilters, filtIDX] = getFilterIDX(filters);
       %% Initialization of motion magnified luma component
    magnifiedLumaFFT = zeros(h,w,nF,'single');
    
    buildLevel = @(im_dft, k) ifft2(ifftshift(croppedFilters{k}.* ...
        im_dft(filtIDX{k,1}, filtIDX{k,2})));
    
    reconLevel = @(im_dft, k) 2*(croppedFilters{k}.*fftshift(fft2(im_dft)));
%%
 %% First compute phase differences from reference frame
    numLevels = numel(filters);        
    fprintf('Moving video to Fourier domain\n');
    vidFFT = zeros(h,w,nF,'single');
    for k = 1:nF
        originalFrame = rgb2ntsc(im2single(vid(:,:,:,k)));
        tVid = imresize(originalFrame(:,:,1), [h w]);
        vidFFT(:,:,k) = single(fftshift(fft2(tVid)));%FFT2 Two-dimensional discrete Fourier Transform.
                                                     %FFTSHIFT Shift zero-frequency component to center of spectrum.           
    end
    clear vid;

    for level = 2:numLevels-1
%     for level = 2
        %% Compute phases of level
        % We assume that the video is mostly static
        pyrRef = buildLevel(vidFFT(:,:,refFrame), level);        
        pyrRefPhaseOrig = pyrRef./abs(pyrRef);
        pyrRef = angle(pyrRef);        

        delta = zeros(size(pyrRef,1), size(pyrRef,2) ,nF,'single');
        fprintf('Processing level %d of %d\n', level, numLevels);
           
        
        for frameIDX = 1:nF
            filterResponse = buildLevel(vidFFT(:,:,frameIDX), level);
            pyrCurrent = angle(filterResponse);
            delta(:,:,frameIDX) = single(mod(pi+pyrCurrent-pyrRef,2*pi)-pi);                          
        end
        
        
%         %% Temporal Filtering
%         fprintf('Bandpassing phases\n');
%         delta = temporalFilter(delta, fl/fs,fh/fs); 


%         %% Apply magnification
% 
%         fprintf('Applying magnification\n');
%         for frameIDX = 1:nF
% 
%             phaseOfFrame = delta(:,:,frameIDX);
%             originalLevel = buildLevel(vidFFT(:,:,frameIDX),level);
%             %% Amplitude Weighted Blur        
%             if (sigma~= 0)
%                 phaseOfFrame = AmplitudeWeightedBlur(phaseOfFrame, abs(originalLevel)+eps, sigma);        
%             end
% 
%             % Increase phase variation
%             phaseOfFrame = magPhase *phaseOfFrame;  
%             
%             if (attenuateOtherFreq)
%                 tempOrig = abs(originalLevel).*pyrRefPhaseOrig;
%             else
%                 tempOrig = originalLevel;
%             end
%             tempTransformOut = exp(1i*phaseOfFrame).*tempOrig; 
% 
%             curLevelFrame = reconLevel(tempTransformOut, level);
%             magnifiedLumaFFT(filtIDX{level,1}, filtIDX{level,2},frameIDX) = curLevelFrame + magnifiedLumaFFT(filtIDX{level,1}, filtIDX{level,2},frameIDX);
%         end
% 
 
    end
    
%% Frequency Spectr of each pixel
deltafreq=delta/(2*pi);


for i=1:size(deltafreq,3)
    yfreq(:,i)=reshape(deltafreq(:,:,i),size(deltafreq,1)*size(deltafreq,2),1);
end

xfreq=linspace(0,vr.Duration,vr.NumberOfFrame);

%%
figure(1);
hold on
for i=1:size(yfreq,1)
plot(xfreq,yfreq(i,:))
end
title('Noisy time domain signal of each pixel')
xlabel('Time')
ylabel('Amplitude')
hold off

figure(2);
hold on 
t=xfreq;
fs=t(end)/(length(t)-1);
srate=1/fs;
ts=t;
for i=1:size(yfreq,1)
z=yfreq(i,:);
Z = fft(z,length(ts));
Pzz = Z.*conj(Z)/length(ts);
f = srate/length(ts)*(0:round(length(ts)/2));
plot(f,Pzz(1:(round(length(ts)/2)+1)))
end
title('Power spectral density of each pixel')
xlabel('Frequency (Hz)')
hold off
%% Frequency Spectr of averaged frame freq
figure(3);
deltafreq=delta/(2*pi);
for i=1:size(deltafreq,3)
    yfreq(i)=mean(mean(deltafreq(:,:,i)));
end

xfreq=linspace(0,vr.Duration,vr.NumberOfFrame);

t=xfreq;
fs=t(end)/(length(t)-1);
srate=1/fs;
ts=t;
z=yfreq;

subplot(2,1,1);
plot(ts,z)
title('Noisy time domain signal of averaged frame freq')
xlabel('Time')
ylabel('Amplitude')

Z = fft(z,length(ts));

Pzz = Z.*conj(Z)/length(ts);
f = srate/length(ts)*(0:round(length(ts)/2));
subplot(2,1,2);
plot(f,Pzz(1:(round(length(ts)/2)+1)))

title('Power spectral density of averaged frame freq')
xlabel('Frequency (Hz)')
%% Sum all frequency spectra
deltafreq=delta/(2*pi);


for i=1:size(deltafreq,3)
    yfreq(:,i)=reshape(deltafreq(:,:,i),size(deltafreq,1)*size(deltafreq,2),1);
end

xfreq=linspace(0,vr.Duration,vr.NumberOfFrame);


figure(4);
subplot(2,1,1);
plot(xfreq,sum(yfreq))
title('Sum of Noisy time domain signals of each pixel')
xlabel('Time')
ylabel('Amplitude')


subplot(2,1,2);
t=xfreq;
fs=t(end)/(length(t)-1);
srate=1/fs;
ts=t;
Pzz=0;

for i=1:size(yfreq,1)
z=yfreq(i,:);
Z = fft(z,length(ts));
Pzz = Z.*conj(Z)/length(ts)+Pzz;
f = srate/length(ts)*(0:round(length(ts)/2));
end

plot(f,Pzz(1:(round(length(ts)/2)+1)))
title('Sum of Power spectral density of each pixel')
xlabel('Frequency (Hz)')



%%
clear vidObj
clear vr
delete 'testvideo.avi'