function idx_ch = charging(im)
% Remove noise
im = medfilt2(im);

% Compute histogram
[counts,binLocations] = imhist(im);
% Cumulative distribution 
Bcsn = cumsum(counts)/sum(counts);        % Normalised Frequency Vector
p=0.9;                                   % Percentile parameter
val = find(Bcsn <= p, 1, 'last');         % gray value where <= 0.9

% % Half-domain conversion to [0 1]
% domain = (binLocations(129:end)-0.5)/0.5;
% % weights from percentile till end taken from half-interval
 counts(1:val)=0;
% W= counts(129:end);
% W= W/sum(W);
% % metric 
% idx_ch = 1-domain'*W;


% New approach
pos = 129;
domain_new = binLocations;
W_new= counts;
W_new(1:pos)=0;
W_new= W_new/sum(W_new);
W_new(isnan(W_new))=0;

% metric 
idx_ch_new = domain_new'*W_new;

% Rescaling
if(idx_ch_new<=0.5)
idx_ch = 1.0;
else
idx_ch = 1- interp1([0.5,1],[0,1],idx_ch_new);
end



