%% Post-processing Stripe Removal Results
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 20.04.2024
%-------------------------------------------------------------------------
% This script evaluates the performance of stripe removal by application of
% different metrics. Calculation is automated on all results found in
% subfolders of the name 'GSR', 'VSNR', 'MDSR' 'WFF'.
%
% The evaluation is stored in a csv-file containing all necessary
% information. Calculations include several metrics proposed in (1).
%
% Comment: Please delete existing csv-files created by this script before
% rerunning. Otherwise, entries are simply appended to the existing file.
%
% (1) Roldan et al., Image quality evaluation for FIB-SEM images (2023)

%% Automatic Accumulation of Results
% Specify a folder and run the process for all result images available!
clearvars dir; addpath('Algorithms')
folder    = uigetdir('Select a "Result" folder.');
system_separator = erase(fullfile(' ',' '),' ');
subfolder = erase(folder,[system_separator,'Result']);

methods = {'GSR','VSNR','MDSR','WFF'};

for i = 1:length(methods)
    current_dir = fullfile(folder,methods{i});
    images = {dir(fullfile(current_dir,'*.tif')).name}; % available results
    if isempty(images)
        break
    end
    
    for j = 1:length(images)
        t   = imfinfo(fullfile(current_dir,images{j})); n = numel(t);
        res = zeros(t(1).Height,t(1).Width,n,'single');
        img = res;
        ideal = res;

        image_name = [regexprep(images{j}, '.*_(.*?)_.*', '$1'),'_',regexprep(images{j}, '.*_(.*)', '$1')];
        ideal_name = insertBefore(image_name,'_prep','_Ideal');
        % Check if ideal exists
        isideal = isfile(fullfile(subfolder,'Preprocessed',ideal_name));
       
        % Read images
        for k=1:n
            res(:,:,k) = single(imread(fullfile(folder,methods{i},images{j}),'Index',k));
            img(:,:,k) = single(imread(fullfile(subfolder,'Preprocessed',image_name),'Index',k));
            if isideal
                ideal(:,:,k) = single(imread(fullfile(subfolder,'Preprocessed',ideal_name),'Index',k));
            end
        end
        Slicestart = str2double(regexprep(image_name,'.*_prep(.*\d+?)-\d+.tif','$1'));
        Sliceend   = str2double(regexprep(image_name,'.*_prep\d+-(.*\d?)+.tif','$1'));

        % Calculation of metrics
        % Comparative metrics (img <-> ideal, result <-> ideal)
        PSNR = NaN(n,2);
        SSIM = NaN(n,2);
        MSSSIM = NaN(n,2);
        RMSE = NaN(n,2);
        
        % Individual metrics (img,ideal,result)
        Curtaining = NaN(n,3);
        Blur = NaN(n,3);
        Charging = NaN(n,3);
        Contrast = NaN(n,3);
        Noise = NaN(n,3);
        
        for k = 1:n
            if isideal
                PSNR(k,:)       = [psnr(img(:,:,k),ideal(:,:,k)),psnr(res(:,:,k),ideal(:,:,k))];
                SSIM(k,:)       = [ssim(img(:,:,k),ideal(:,:,k)),ssim(res(:,:,k),ideal(:,:,k))];
                MSSSIM(k,:)     = [multissim(img(:,:,k),ideal(:,:,k)),multissim(res(:,:,k),ideal(:,:,k))];
                RMSE(k,:)       = [sqrt(mean((img(:,:,k) - ideal(:,:,k)).^2,'all')),sqrt(mean((res(:,:,k) - ideal(:,:,k)).^2,'all'))];
                Curtaining(k,:) = [band_curtaining(img(:,:,k)),band_curtaining(ideal(:,:,k)),band_curtaining(res(:,:,k))];
                Blur(k,:)       = [BlurLevel(img(:,:,k)),BlurLevel(ideal(:,:,k)),BlurLevel(res(:,:,k))];
                Charging(k,:)   = [charging(img(:,:,k)),charging(ideal(:,:,k)),charging(res(:,:,k))];
                Contrast(k,:)   = [contrastMetric(img(:,:,k)),contrastMetric(ideal(:,:,k)),contrastMetric(res(:,:,k))];
                Noise(k,:)      = [NoiseLevel(img(:,:,k)),NoiseLevel(ideal(:,:,k)),NoiseLevel(res(:,:,k))];
            else
                Curtaining(k,:) = [band_curtaining(img(:,:,k)),NaN,band_curtaining(res(:,:,k))];
                Blur(k,:)       = [BlurLevel(img(:,:,k)),NaN,BlurLevel(res(:,:,k))];
                Charging(k,:)   = [charging(img(:,:,k)),NaN,charging(res(:,:,k))];
                Contrast(k,:)   = [contrastMetric(img(:,:,k)),NaN,contrastMetric(res(:,:,k))];
                Noise(k,:)      = [abs(NoiseLevel(img(:,:,k))),NaN,abs(NoiseLevel(res(:,:,k)))];
            end
        end

	    % Create and write table
        convertedname =  regexprep(images{j}, '(\d)p(\d)', '$1.$2');
        numbers = regexp(convertedname, '\d+(\.\d+)?', 'match');
        % Set string of parameters
        parameters = '';
        if strcmp(methods{i},'GSR')
            if str2double(numbers{1}) == 2
                method = 'GSR 2D';
                parameters = ['mu=(',numbers{2},',',numbers{3},'), iterations=',numbers{4},', proj=',numbers{5},', normalize=',numbers{6}];
            else
                method = 'GSR 3D';
                parameters = ['mu=(',numbers{2},',',numbers{3},'), iterations=',numbers{4},', proj=',numbers{5},', resz=',numbers{6},', normalize=',numbers{7}];
            end
        elseif strcmp(methods{i},'MDSR')
            method = 'MDSR';
            parameters = ['dir=',numbers{1},', dec=',numbers{2},', sigma=',numbers{3},', sigma_a=',numbers{4}];
        elseif strcmp(methods{i},'VSNR')
            method = 'VSNR';
            parameters = ['alpha=(',numbers{1},',',numbers{2},',',numbers{3},'), steps=',numbers{4}];
        elseif strcmp(methods{i},'WFF')
            method = 'WFF';
            parameters = ['dec=',numbers{1},', wavelet=db',numbers{2},', sigma=',numbers{3}];
        end

        % Create data table
        varnames = {'Image name',...
                    'Slice',...
                    'Method',...
                    'Parameters',...
                    'PSNR(Image,Ideal)','PSNR(Result,Ideal)',...
                    'SSIM(Image,Ideal)','SSIM(Result,Ideal)',...
                    'MS-SSIM(Image,Ideal)','MS-SSIM(Result,Ideal)',...
                    'RMSE(Image,Ideal)','RMSE(Result,Ideal)',...
                    'Curtaining(Image)','Curtaining(Ideal)','Curtaining(Result)',...
                    'Blur(Image)','Blur(Ideal)','Blur(Result)',...
                    'Charging(Image)','Charging(Ideal)','Charging(Result)',...
                    'Contrast(Image)','Contrast(Ideal)','Contrast(Result)',...
                    'Noise Level(Image)','Noise Level(Ideal)','Noise Level(Result)'};
        T = table(repmat({extractBefore(image_name,'_prep')},n,1),...
                  (Slicestart:Sliceend)', ...
                  repmat({method},n,1),...
                  repmat({parameters},n,1),...
                  PSNR(:,1),PSNR(:,2),...
                  SSIM(:,1),SSIM(:,2),...
                  MSSSIM(:,1),MSSSIM(:,2),...
                  RMSE(:,1),RMSE(:,2),...
                  Curtaining(:,1),Curtaining(:,2),Curtaining(:,3),...
                  Blur(:,1),Blur(:,2),Blur(:,3),...
                  Charging(:,1),Charging(:,2),Charging(:,3),...
                  Contrast(:,1),Contrast(:,2),Contrast(:,3),...
                  Noise(:,1),Noise(:,2),Noise(:,3),'VariableNames',varnames);

        if isfile(fullfile(folder,[erase(image_name,'.tif'),'-Evaluation.csv']))
            writetable(T ,fullfile(folder,[erase(image_name,'.tif'),'-Evaluation.csv']),'WriteMode','append','Delimiter',';')
        else
            writetable(T ,fullfile(folder,[erase(image_name,'.tif'),'-Evaluation.csv']),'Delimiter',';')
        end
    end
end