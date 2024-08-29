%% Helper Function: LoadImage
%-------------------------------------------------------------------------
% Author: Niklas Rottmayer
% Date:   18.06.2024
%-------------------------------------------------------------------------
%
% Loads an image (2D/3D) automatically
% Input: 
% - path - full path to the image file
% Output: 
% - img  - 2D/3D image
% - n    - depth of the image stack
%

function [img,n] = LoadImage(path)
    % Read image data
    t_image = imfinfo(path); n = numel(t_image);
    if n == 1
        img = single(imread(path));
    elseif n > 1
        img = zeros(t_image(1).Height,t_image(1).Width,n,'single');
        for i=1:n
            img(:,:,i) = single(imread(path,'Index',i));
        end
    end
end