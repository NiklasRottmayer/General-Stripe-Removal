%% Preprocessing for Stripe Removal
%
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date: 17.04.2024
%-------------------------------------------------------------------------
%
% This script performs pre-processing on an image (2D or 3D) to prepare for
% application of stripe removal methods with the following steps:
%
% 1. Normalization      -   transform value range into [0,1]    
%    (required)             
%
% 2. Permutation        -   permute image dimension manunally such that 
%    (optional)             stripes are vertical (in y/2nd-direction)
%
% 3. Slice-extraction   -   select slices of a 3D image to process 
%    (optional)
%
% 4. Storage            -   save pre-processed image accordingly
%    (required)
% 
% Comment: It is recommended to perform normalization first, i.e., on the
% entire image stack. The name of an 'ideal' image should be '[image
% name]_Ideal.[file format]' for automatic loading. 

%% Step 0: Loading image data (required)
clearvars img ideal

% Select an image from a directory
[image_name,folder] = uigetfile({'*.png; *.jpg; *.jpeg; *.bmp; *.tif; *.tiff','Image Files (*.png,*.jpg,*.jpeg,*.bmp,*.tif,*.tiff)'});

% Read image data
t_image = imfinfo([folder,image_name]); n = numel(t_image); Slices = 1:n;
if n == 1
    img = single(imread([folder,image_name]));
elseif n > 1
    img = zeros(t_image(1).Height,t_image(1).Width,n,'single');
    for i=1:n
        img(:,:,i) = single(imread([folder,image_name],'Index',i));
    end
end

% If available load an 'ideal' image (autmatically)
[~,~,file_type] = fileparts(image_name);
ideal_name = strrep(image_name,file_type,['_Ideal',file_type]);
if isfile([folder,ideal_name])
    t_ideal = imfinfo([folder,ideal_name]); n_ideal = numel(t_image);
    if n_ideal == 1
        ideal = single(imread([folder,ideal_name]));
    elseif n_ideal > 1 
        ideal = zeros(t_ideal(1).Height,t_ideal(1).Width,n_ideal,'single');
        for i=1:n
            ideal(:,:,i) = single(imread([folder,image_name],'Index',i));
        end
    end
end

if contains(image_name,'_')
    image_name = strrep(image_name,'_','-');
    ideal_name = strrep(image_name,file_type,['_Ideal',file_type]); % update
    warning(['The image name has been modified by replacing "_" with "-". ' ...
             'This is only necessary for a correct functionality when post-processing.'])    
end

%% Step 1: Normalization to [0,1] (required)
% Applies normalization to values [0,1] for img
% If available, the ideal is normalized with the same values, i.e. its
% values are not necessarily in [0,1]!
immax = max(img,[],'all');
immin = min(img,[],'all');
img = (img - immin)/(immax - immin);
if exist('ideal','var')
    ideal = (ideal - immin)/(immax - immin);
end

%% Step 2: Permutate dimensions (optional)
% Permute dimensions such that stripes are vertical (2nd/y-direction)
dir = [2,1,3];
img = permute(img,dir);
if exist('ideal','var')
    ideal = permute(ideal,dir);
end

%% Step 3: Slice extraction (optional)
% Extract selected slices or bands of slices
Slices = 60:60; 
if max(Slices) > n || min(Slices) < 1
    Slices = 1:n;
    error('Please make a valid selection of slices.')
end

img = img(:,:,Slices);
if exist('ideal','var')
    ideal = ideal(:,:,Slices);
end

%% Step 4: Data storage (required)
% Save the pre-processed images in a subsequent folder called
% "Preprocessed" placed in the same directory as the image. 
% Images are stored as tif-images to retain the normalized values.

if ~isfolder([folder,'Preprocessed'])
    mkdir([folder,'Preprocessed'])
end

if isequal(Slices,1:n)
    save_name       = fullfile(folder,'Preprocessed',regexprep(image_name, '\.(png|jpg|jpeg|bmp|tif|tiff)','_prep.tif'));
    save_ideal_name = fullfile(folder,'Preprocessed',regexprep(ideal_name, '\.(png|jpg|jpeg|bmp|tif|tiff)','_prep.tif'));
elseif ~isequal(Slices,1:n)
    save_name       = fullfile(folder,'Preprocessed',regexprep(image_name, '\.(png|jpg|jpeg|bmp|tif|tiff)',['_prep',num2str(min(Slices,[],'all')),'-',num2str(max(Slices,[],'all')),'.tif']));
    save_ideal_name = fullfile(folder,'Preprocessed',regexprep(ideal_name, '\.(png|jpg|jpeg|bmp|tif|tiff)',['_prep',num2str(min(Slices,[],'all')),'-',num2str(max(Slices,[],'all')),'.tif']));
end

SaveImage(img,save_name);
if exist('ideal','var')
    SaveImage(ideal,save_ideal_name)
end

