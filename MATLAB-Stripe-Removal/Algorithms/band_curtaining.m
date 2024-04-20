%% Curtaining Metric 
% Function to calculate the curtaining metrics proposed by Roldan:
% [Roldan et al., Image quality evaluation for FIB-SEM images (2023)]

function index_c = band_curtaining(im)
im = im(:,:,1);
% 1. crop the image to be square 
sz = size(im);
im = im(1:min(sz),1:min(sz));
sz = size(im);
mImage = sz(1,1);
nImage = sz(1,2);

[Gx, Gy] = imgradientxy(im);

% 1. compute the FFT
fGx = fft2(Gx);
fdx = fftshift(fGx);
% Visualize
fdx = log(1 + abs(fdx));
fdx = fdx/max(max(fdx));
 

%3.1 select band analysis
% Crop fdx in x and y

% New-region selection 121x31
x_right =ceil(nImage/2)+60;
x_left  =ceil(nImage/2)-60;
y_up =ceil(nImage/2)+15;
y_bottom  =ceil(nImage/2)-15;


fdx = fdx(y_bottom:y_up,x_left:x_right);

sz_fdx = size(fdx);
fd_binary = zeros(sz_fdx(1),sz_fdx(2));

% Detect maxima per column
for i=1:sz_fdx(2)
[val pos] = max(fdx(:,i));
fd_binary(pos,i) =1;
end
signal = sum(fd_binary')/sz_fdx(2);
[v_center posit] = max(signal);
if((posit==1)||(posit==length(signal)))
    posit = round(length(signal)*0.5);
end
index_c_aux = 1-sum(signal(posit-1:posit+1));

% Rescaling
if(index_c_aux>0.9)
index_c = 1.0;
else
index_c = interp1([0,0.9],[0,1],index_c_aux);
end

end