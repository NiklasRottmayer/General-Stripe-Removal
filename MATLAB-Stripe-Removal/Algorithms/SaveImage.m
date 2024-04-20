%% Helper Function: Save Image

function SaveImage(img,filename,BitsPerSample)
    if nargin < 3
        BitsPerSample = 32;
    end

    if contains(filename,{'.png','.jpg','.jpeg','.bmp'})
        imwrite(img,filename)
    else
        t = Tiff(filename,'w');
        %
        tagstruct.ImageLength = size(img, 1);
        tagstruct.ImageWidth = size(img, 2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.Photometric = Tiff.Photometric.LinearRaw;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.BitsPerSample = BitsPerSample; % 32;
        tagstruct.SamplesPerPixel = 1; % Color depth - grayscale
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        %
        t.setTag(tagstruct);
        if size(img,3) > 1
            for l = 1:size(img,3)
                t.write(squeeze(img(:,:,l)));
                if l<size(img,3)
                    writeDirectory(t);
                    t.setTag(tagstruct);
                end
            end
        else 
            t.write(img);
        end
        t.close();
    end
end
