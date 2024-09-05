%% Processing Script for Stripe Removal
%
% This script performs stripe removal for corrupted images.
% Ideally, images should be pre-processed using the 'PreprocessingScript'.
% The following methods are supported:
% 1. (ours) GeneralStripeRemover: A variational method that penalizes 
%    undesired image and stripe features to remove artifacts. It supports
%    processing in 3D and extends on the work by Liu et al. [1].
% 2. MDSR*: A powerful Fourier filtering approach based on the
%    non-subsampled contourlet transform, a multi-scale and directional
%    decomposition. Only slice-wise processing is currently supported.
%    Compared to the original proposition [2] we restrict filtering to
%    subimages of a limited range of directions.
% 3. VSNR: A versatile variational method capable of removing repeating
%    artifacts through the choice of patterns. Only slice-wise processing
%    is currently supported [3].
% 4. WFF: Wavelet Fourier-filtering is one of the earliest successfull
%    implementations of Fourier filtering for removing stripes based on the
%    wavelet decomposition. It is commonly referenced and compared with in
%    literature. [4]
%
% [1] Liu et al., Oblique stripe removal in remote sensing images via
%     oriented variation (2018).
% [2] Liang et al., Stripe artifact elimination based on nonsubsampled 
%     contourlet transform for light sheet fluorescence microscopy (2016).
% [3] Fehrenbach et al., Variational algorithms to remove stationary noise: 
%     Applications to microscopy imaging (2012).
% [4] Münch et al., Stripe and ring artifact removal with combined wavelet—
%     Fourier filtering (2009).
%
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 29.08.2024
%-------------------------------------------------------------------------
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
if n == 1, img = single(imread([folder,image_name]));
elseif n > 1
    img = zeros(t_image(1).Height,t_image(1).Width,n,'single');
    for i=1:n, img(:,:,i) = single(imread([folder,image_name],'Index',i)); end
end

system_separator = erase(fullfile(' ',' '),' ');

% Create storage folder for results
if contains(folder,[system_separator,'Preprocessed',system_separator])
    result_folder = strrep(folder,[system_separator,'Preprocessed',system_separator],[system_separator,'Result',system_separator]);
else
    result_folder = fullfile(folder,'Result',system_separator);
end
if ~isfolder(result_folder), mkdir(result_folder), end

%% Method 1: (ours) General Stripe Remover
% Interpretation
% - mu1 -> strength of stripe removal
% - mu2 -> precaution of removing structures
%
% Suggested parameters mu:
% - [5, 0.1]  or [7, 0.1] if stripes are thin and impairment is low.  
% - [10, 0.1] or [12,0.2] if stripes are wider and corruptions severely 
%                         influence the visual impression.
% - [15, 0.5]             if corruptions are severe and stripes are of
%                         short length (on the scale of structures).
% In some instances a more conservative removal using [3,0.05] also led to
% great results. The suggested settings might not yield "optimal" results
% but show the range of values. Compared to the paper mu is rescaled by a
% factor of 30. 

% Parameters
mu = [7,0.1];           % Weighting parameters [mu1,mu2] (default: [10,0.1])
steps = 25000;          % Number of iterations (default: 25000, >=5000 suggested)
resz = 1;               % Resolution in z-direction (default: 1)
proj = 1;               % 1 for GSR method in the paper, 0 for [Liu et. al 2018]
normalize = 0;          % 1 to normalize input, 0 for direct processing (default: 0)
direction = [1,0,0];    % Stripe direction (array[3])
verbose = 1;            % 1 for step counter, 0 for no outputs

if ~isfolder(fullfile(result_folder,'GSR')), mkdir(fullfile(result_folder,'GSR')), end

assert(norm(direction,1) > 0), dimension = 3; direction = direction./norm(direction,1);
if ismatrix(img) || resz == 0, resz = 0; dimension = 2; direction = direction(1:2)./norm(direction(1:2),1); end
% Processing
u = GeneralStripeRemover(img,steps,mu,proj,resz,normalize,direction,verbose);
% Storage
save_name = fullfile(result_folder,'GSR',['GSR',num2str(dimension),'D_mu',strrep(num2str(mu(1)),'.','p'),'-',...
                    strrep(num2str(mu(2)),'.','p'),'_steps',strrep(num2str(steps),'.','p'),...
                    '_proj',num2str(proj),'_resz',strrep(num2str(resz),'.','p'),...
                    '_direction',strjoin(arrayfun(@(x) regexprep(regexprep(strrep(sprintf('%.2f',round(x, 2)),'.','p'),'p?0*$','p'),'p$',''),direction,'UniformOutput',false), '-'),...
                    '_normalize',num2str(normalize),'_',image_name]);
SaveImage(u,save_name);

%% Method 2: Multidirectional Stripe Remover (MDSR*) [2]
% Suggested parameters: 
% dirnum = 8, decnum = 5, sigma in [5,25], sigmaa = 0.3, 
% max_angle = pi/4
% Comment: Only sigma needs to be tuned. decnum < 5 only when image is too 
% small and an error is raised. 

% Parameters
sigma = 10;                 % Gaussian damping parameter (default: 10) 
dirnum = 8;                 % Number of directions (default: 8)
decnum = 5;                 % Decomposition depths (default: 5)            
sigmaa = 0.3;               % Damping fall-off parameter (default: 0.3)
max_angle = pi/4;           % Restrict filtering to directions pi/2+-max_angle (default: pi/4)
                            % where pi/2 is the vertical direction.
                            % max_angle >= pi/2 will result in the original method.

if ~isfolder(fullfile(result_folder,'MDSR')), mkdir(fullfile(result_folder,'MDSR')), end

% Processing
u = zeros(size(img),'single');
for i = 1:n
    u(:,:,i) = single(MDSR(double(img(:,:,i)),dirnum,decnum,sigma,sigmaa,max_angle));
end 
% Storage
save_name = fullfile(result_folder,'MDSR',['MDSR_dir',strrep(num2str(dirnum),'.','p'),...
                '_dec',strrep(num2str(decnum),'.','p'),...
                '_sigma',strrep(num2str(sigma),'.','p'),'_sigmaa',strrep(num2str(sigmaa),'.','p'),'_',strrep(num2str(sigma),'.','p'),'_angle',strrep(num2str(round(max_angle,2)),'.','p'),'_',image_name]);
SaveImage(u,save_name);

%% Method 3: Variational Stationary Noise Remover (VSNR) [3]
% Suggested parameters:
% alpha in [1,13]^3 depending on appearance of thin, medium and large
% stripes. 
% p = 1 and steps = 25000

% 2D Artifact patterns
Filters = [];
Filters(:,:,1) = GaborFilter(size(img),0.1,0,0,0.5,0.1); % Small 
Filters(:,:,2) = GaborFilter(size(img),0.1,0,0,1,0.1); % Medium
Filters(:,:,3) = GaborFilter(size(img),0.1,0,0,2,0.1); % Large

% Parameters
alpha = [5,5,5];    % weighting parameters (default: [5,5,5])           
p = 1;              % value or array specifying the p-norm used for penalization (default: 1)
steps = 25000;       % Number of iterations (default: 25000)

% Further parameters - for details see [3]
eps = 10e-3; prec = 10e-3; C = 1; % no change necessary

if ~isfolder(fullfile(result_folder,'VSNR')), mkdir(fullfile(result_folder,'VSNR')), end

% Processing
u = zeros(size(img),'single');
for i = 1:n
    u(:,:,i) = VSNR(img(:,:,i),eps,p,Filters,alpha,steps,prec,C);
end
%Storage
save_name = fullfile(result_folder,'VSNR','VSNR_alpha');
for i = 1:size(Filters,3)-1, save_name = [save_name,strrep(num2str(alpha(i)),'.','p'),'-']; end
save_name = [save_name,strrep(num2str(alpha(end)),'.','p'),'_steps',strrep(num2str(steps),'.','p'),'_',image_name];
SaveImage(u,save_name);

%% Method 4: Wavelet Fourier filtering (WFF) [4]
% Suggested parameters:
% decnum = 6, wname = 'db20' and sigma in [5,25]
% Comment: decnum < 6 only for small images when error is raised.

% Parameters
decnum = 6;    % Decomposition depths (default: 6)
wname = 'db20'; % Wavelet type (default: 'db20')
sigma = 10;     % Gaussian damping parameter (default: 10)

if ~isfolder(fullfile(result_folder,'WFF')), mkdir(fullfile(result_folder,'WFF')), end

% Processing
u = zeros(size(img),'single');
for i = 1:n
    u(:,:,i) = WFF(img(:,:,i),decnum,wname,sigma);
end
% Storage
save_name = fullfile(result_folder,'WFF',['WFF_dec',strrep(num2str(decnum),'.','p'),...
                '_wname-',wname,'_sigma',strrep(num2str(sigma),'.','p'),'_',image_name]);
SaveImage(u,save_name);
