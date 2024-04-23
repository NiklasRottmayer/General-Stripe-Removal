%% Processing Stripe Removal
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 17.04.2024
%-------------------------------------------------------------------------
% 
% This script performs stripe removal for corrupted images.
% Ideally, images should be pre-processed using the 'PreprocessingScript'.
% The following methods are supported:
%
% 1. GSR: The general stripe remover is an extension of the variational
% method proposed by Liu et al. in 2018. It utilizes basic penalizations
% terms for stripe removal combined with L1-regularization on the stripe
% image. It is the only method which supports non-slice wise 3D processing.
% [Liu et al., Oblique stripe removal in remote sensing images via
%  oriented variation (2018)]
%
% 2. MDSR: A powerful Fourier filtering approach based on the
% non-subsampled contourlet transform for a multi-scale and -directional
% decomposition. Images are processed slice-wise.
% NOTICE: Requires the 'Nonsubsampled Contourlet Transform' (Add-on) and 
%         needs double inputs. Read MDSR.m installation help.
% [Liang et al., Stripe artifact elimination based on nonsubsampled 
%  contourlet transform for light sheet fluorescence microscopy (2016)]
%
% 3. VSNR: A versatile variational method which can be adapted for general
% repetitive artifact patterns. It heavily relies on the decomposition of
% artifacts as convolution of locations and patterns. The methods performs
% 2D processing, i.e. slice-wise processing for 3D.
% [Fehrenbach et al., Variational algorithms to remove stationary noise: 
%  Applications to microscopy imaging (2012)]
%
% 4. WFF: Wavelet Fourier-filtering is one of the earliest successfull
% implementations of Fourier filtering for removing stripes. It relies on
% the wavelet decomposition for a multi-directional and -scale separation
% of the image. 
% [Münch et al., Stripe and ring artifact removal with combined wavelet —
%  Fourier filtering (2009)]
%
%% Step 0: Loading image data (2D or 3D)
addpath('Algorithms')

% Select an image from your directory
[image_name,folder] = uigetfile({'*.png; *.jpg; *.jpeg; *.tif; *.tiff','Image Files (*.png,*.jpg,*.jpeg,*.tif,*.tiff)'});

% Check if pre-processing has been performed
if ~contains(image_name,'_prep')
    warning(['The selected image has not been pre-processed using the "PreprocessingScript.m".\n' ...
             'Please make sure, that\n' ...
             '1. the image is normalized to [0,1]\n' ...
             '2. stripe are vertical.\n'])
end

% Load image data
t_image = imfinfo([folder,image_name]); n = numel(t_image);
if n == 1
    img = single(imread([folder,image_name]));
elseif n > 1
    img = zeros(t_image(1).Height,t_image(1).Width,n,'single');
    for i=1:n
        img(:,:,i) = single(imread([folder,image_name],'Index',i));
    end
end

system_separator = erase(fullfile(' ',' '),' ');

% Create storage folder for results
if contains(folder,[system_separator,'Preprocessed',system_separator])
    result_folder = strrep(folder,[system_separator,'Preprocessed',system_separator],[system_separator,'Result',system_separator]);
else
    result_folder = fullfile(folder,'Result',system_separator);
end
if ~isfolder(result_folder)
    mkdir(result_folder)
end

%% Method 1: General Stripe Remover [Rottmayer et al. 2024, Liu et. al 2018]
% Weighting parameters suggestions: [10,0.1], [15,0.5], [12,0.2], [7,0.1]
% see "GeneralStripeRemover3D.m" or "GeneralStripeRemover2D.m". Compared to
% the paper parameters are scaled by factor 30.

% Parameters
mu = [10,0.1];  % Weighting parameters [mu1,mu2] (default: [10,0.1])
steps = 25000;  % Number of iterations (default: 25000)
resz = 1;       % Resolution in z-direction (if applicable) (default: 1)
proj = 1;       % 1 for GSR method in the paper, 0 for [Liu et. al 2018]
normalize = 0;  % 1 to normalize input, 0 for direct processing (default: 0)
GPU = 1;        % 1 for use of GPU (if available), 0 for use of CPU
verbose = 1;    % 1 for step counter, 0 for no outputs
D3 = 1;         % 1 for processing 3D, 0 for slice-wise processing (only applicable for 3D images)

if ~isfolder(fullfile(result_folder,'GSR'))
    mkdir(fullfile(result_folder,'GSR'))
end
if D3 && n > 1
    u = GSR3D(img,steps,mu,proj,resz,normalize,GPU,verbose);
    save_name = fullfile(result_folder,'GSR',['GSR3D_mu',strrep(num2str(mu(1)),'.','p'),'-',...
                    strrep(num2str(mu(2)),'.','p'),'_steps',strrep(num2str(steps),'.','p'),...
                    '_proj',num2str(proj),'_resz',strrep(num2str(resz),'.','p'),...
                    '_normalize',num2str(normalize),'_',image_name]);
    SaveImage(u,save_name);
else
    u = GSR2D(img,steps,mu,proj,normalize,GPU,verbose);
    save_name = fullfile(result_folder,'GSR',['GSR2D_mu',strrep(num2str(mu(1)),'.','p'),'-',...
                    strrep(num2str(mu(2)),'.','p'),'_steps',strrep(num2str(steps),'.','p'),...
                    '_proj',num2str(proj),'_normalize',num2str(normalize),'_',image_name]);
    SaveImage(u,save_name);
end

%% Method 2: Multidirectional Stripe Remover (MDSR) [Liang et. al 2016]*
% Parameters suggestions: dirnum = 8, decnum = 5, sigma in [10,15,20]
% * We restrict filtering to a subset of directional subimages by angle.

% Parameters
dirnum = 8;                 % Number of directions (default: 8)
decnum = 5;                 % Decomposition depths (default: 5)
sigma = 10;                 % Gaussian damping parameter (default: 10)
sigmaa = 0.3;               % Damping fall-off parameter (default: 0.3)
max_angle = pi/4;           % Restrict filtering to directions pi/2+-max_angle (default: pi/4)
                            % where pi/2 is the vertical direction.
                            % max_angle >= pi/2 will result in the original method.

if ~isfolder(fullfile(result_folder,'MDSR'))
    mkdir(fullfile(result_folder,'MDSR'))
end

u = zeros(size(img),'single');
for i = 1:n
    u(:,:,i) = single(MDSR(double(img(:,:,i)),dirnum,decnum,sigma,sigmaa,max_angle));
end 

save_name = fullfile(result_folder,'MDSR',['MDSR_dir',strrep(num2str(dirnum),'.','p'),...
                '_dec',strrep(num2str(decnum),'.','p'),...
                '_sigma',strrep(num2str(sigma),'.','p'),'_sigmaa',strrep(num2str(sigmaa),'.','p'),'_',image_name]);
SaveImage(u,save_name);

%% Method 3: Variational Stationary Noise Remover (VSNR) [Fehrenbach et al. 2012]
% 2D Artifact patterns
Filters = [];
Filters(:,:,1) = GaborFilter(size(img),0.1,0,0,0.5,0.1); % Small 
Filters(:,:,2) = GaborFilter(size(img),0.1,0,0,1,0.1); % Medium
Filters(:,:,3) = GaborFilter(size(img),0.1,0,0,2,0.1); % Large

% Parameters
alpha = [5,5,5];    % weighting parameters (default: [5,5,5])           
p = 1;              % value or array specifying the p-norm used for penalization (default: 1)
steps = 25000;      % Number of iterations (default: 25000)
GPU = 0;            % 1 for use of GPU (if available), 0 for use of CPU

% Further parameters - for details see [Fehrenbach et al. 2012]
eps = 10e-3; prec = 10e-3; C = 1; % no change necessary


u = zeros(size(img),'single');
for i = 1:n
    u(:,:,i) = fast_VSNR(img(:,:,i),eps,p,Filters,alpha,steps,prec,C,GPU);
end
save_name = fullfile(result_folder,'VSNR','VSNR_alpha');
for i = 1:size(Filters,3)-1
    save_name = [save_name,strrep(num2str(alpha(i)),'.','p'),'-'];
end
if ~isfolder(fullfile(result_folder,'VSNR'))
    mkdir(fullfile(result_folder,'VSNR'))
end
save_name = [save_name,strrep(num2str(alpha(end)),'.','p'),'_steps',strrep(num2str(steps),'.','p'),'_',image_name];
SaveImage(u,save_name);

%% Method 4: Wavelet Fourier filtering (WFF) [Münch et al. 2009]
% Parameters
decnum = 6;    % Decomposition depths (default: 6)
wname = 'db20'; % Wavelet type (default: 'db20')
sigma = 10;     % Gaussian damping parameter (default: 10)

u = zeros(size(img),'single');
for i = 1:n
    u(:,:,i) = WFF(img(:,:,i),decnum,wname,sigma);
end

if ~isfolder(fullfile(result_folder,'WFF'))
    mkdir(fullfile(result_folder,'WFF'))
end
save_name = fullfile(result_folder,'WFF',['WFF_dec',strrep(num2str(decnum),'.','p'),...
                '_wname-',wname,'_sigma',strrep(num2str(sigma),'.','p'),'_',image_name]);
SaveImage(u,save_name);
