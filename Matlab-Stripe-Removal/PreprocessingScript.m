%% Preprocessing for Stripe Removal
%
% This script performs pre-processing on 2D or 3D images to prepare for the
% application of different stripe removal methods. It is highly recommended
% to run this script before proceeding to remove stripes. The following
% steps are supported:
% 1. Normalization      -   rescales image values to the interval [0,1]
%    (recommended)
% 2. Permutation        -   change of image dimensions such that stripes
%    (optional)             are vertical.  
% 3. Slicing            -   selection of slices / substack from 3D image
%    (optional)
% 4. Saving             -   save pre-processing result
%    (recommended)
%
% Comment: If a ground truth to an image exists in the same folder as the
% image, it should be named in the format '[image name]_Ideal.[file format]'
% for automatic loading. 
%
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 29.08.2024
%-------------------------------------------------------------------------

%% Step 0: Loading image data (required)
clearvars img ideal
addpath('Algorithms')

% Select an image from a directory
[image_name,folder] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif; *.tiff','Image Files (*.png,*.jpg,*.jpeg,*.bmp,*.tif,*.tiff)'});

% Read image data
t_image = imfinfo([folder,image_name]); n = numel(t_image); Slices = 1:n;
if n == 1,  img = single(imread([folder,image_name]));
elseif n > 1
    img = zeros(t_image(1).Height,t_image(1).Width,n,'single');
    for i=1:n, img(:,:,i) = single(imread([folder,image_name],'Index',i)); end
end

% If available load an 'ideal' image (autmatically)
[~,~,file_type] = fileparts(image_name);
ideal_name = strrep(image_name,file_type,['_Ideal',file_type]);
if isfile([folder,ideal_name])
    t_ideal = imfinfo([folder,ideal_name]); n_ideal = numel(t_image);
    if n_ideal == 1, ideal = single(imread([folder,ideal_name]));
    elseif n_ideal > 1 
        ideal = zeros(t_ideal(1).Height,t_ideal(1).Width,n_ideal,'single');
        for i=1:n
            ideal(:,:,i) = single(imread([folder,ideal_name],'Index',i));
        end
    end
end

if contains(image_name,'_')
    image_name = strrep(image_name,'_','-');
    ideal_name = strrep(image_name,file_type,['_Ideal',file_type]); % update
    warning(['The image name has been modified by replacing "_" with "-". ' ...
             'This is only necessary for a correct functionality when post-processing.'])    
end

%% Step 1: Normalization to [0,1] (recommended)
% If available the ground truth is normalized by the same values, i.e. not
% necessarily in [0,1] to preserve comparability. 
immax = max(img,[],'all'); immin = min(img,[],'all'); 
img = rescale(img,0,1);
if exist('ideal','var'), ideal = (ideal - immin)/(immax - immin); end

%% Step 2: Permutate dimensions (optional)
% Changes order of dimensions.
dir = [2,1,3]; 
img = permute(img,dir);
if exist('ideal','var'), ideal = permute(ideal,dir); end

%% Step 3: Slice extraction (optional)
% Extract selected slices or bands of slices.
Slices = 1:50; 
Slices = Slices(Slices <= n & Slices >= 1); % only valid selection

img = img(:,:,Slices);
if exist('ideal','var'), ideal = ideal(:,:,Slices); end

%% Step 4: Data storage (required)
% Save the pre-processed images in a subsequent folder called
% "Preprocessed" placed in the same directory as the image. 
% Images are stored as tif-images to retain the normalized values.

if ~isfolder([folder,'Preprocessed']), mkdir([folder,'Preprocessed']), end

if isequal(Slices,1:n)
    save_name       = fullfile(folder,'Preprocessed',regexprep(image_name, '\.(png|jpg|jpeg|bmp|tif|tiff)','_prep.tif'));
    save_ideal_name = fullfile(folder,'Preprocessed',regexprep(ideal_name, '\.(png|jpg|jpeg|bmp|tif|tiff)','_prep.tif'));
elseif ~isequal(Slices,1:n)
    save_name       = fullfile(folder,'Preprocessed',regexprep(image_name, '\.(png|jpg|jpeg|bmp|tif|tiff)',['_prep',num2str(min(Slices,[],'all')),'-',num2str(max(Slices,[],'all')),'.tif']));
    save_ideal_name = fullfile(folder,'Preprocessed',regexprep(ideal_name, '\.(png|jpg|jpeg|bmp|tif|tiff)',['_prep',num2str(min(Slices,[],'all')),'-',num2str(max(Slices,[],'all')),'.tif']));
end

SaveImage(img,save_name);
if exist('ideal','var'), SaveImage(ideal,save_ideal_name), end

