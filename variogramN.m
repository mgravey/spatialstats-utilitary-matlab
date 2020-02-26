function [gh11, nh11] = variogramN(x1)

% function to compute variograms or covariograms, in 1D or 2D,
% the data are on a (possibly incomplete) regular grid.
% the program computes variograms in the frequency domain using
% 2D-FFT-
%
% INPUT:       x1: data matrix. Missing values are indicated by NaN
%           icode: a code to indicate which function to compute
%              =1: variogram
%              =2: covarigram
%
% OUTPUT:    gh11: variogram or covariogram depending on icode
%            nh11: number of pairs available
%
% this program uses the function FFT2, IFF2, FFT2SHIFT and CONJ
% which are standard MATLAB functions.
%
% MORE INFO -> (MARCOTTE, D. 1996) paper in Computers & Geosciences

s=size(x1);

% find the closest multiple of 8 to obtain a good compromise between
% speed (a power of 2) and memoty required
sr2=ceil((2*s-1)/8)*8;

% form a indicator matrix: 1's for all data values
%                          0's for missing values
% in data matrix, replace missing values by 0;

x1id       = ~isnan(x1);                 % 1 for a data values; 0 for missing
x1(~x1id)  = 0;  % missing replaced by 0

fx1        = fftn(x1, sr2);         % Fourier transform of x1


fx1_x1 = fftn(x1 .* x1, sr2);   % Fourier transform of x1*x1


fx1id      = fftn(x1id, sr2);       % Fourier transform of the indicator matrix


% compute number of pairs at all lags
nh11       = round(real(ifftn(conj(fx1id) .* fx1id)));

% compute the different structural functions according to icode

                          % variogram is computed
gh11   = real(ifftn(conj(fx1id) .* fx1_x1 + conj(fx1_x1) .* fx1id-2 * conj(fx1) .* fx1));
gh11   = gh11 ./ max(nh11,1) / 2;



% reduce matrix to required size and shift so that the 0 lag appears at the
% center of each matrix

% nh11 = [nh11(1:n,1:p) nh11(1:n,nc2-p+2:nc2); nh11(nr2-n+2:nr2,1:p) nh11(nr2-n+2:nr2,nc2-p+2:nc2)];
% gh11 = [gh11(1:n,1:p) gh11(1:n,nc2-p+2:nc2); gh11(nr2-n+2:nr2,1:p) gh11(nr2-n+2:nr2,nc2-p+2:nc2)];

gh11 = fftshift(gh11);
nh11 = fftshift(nh11);

dims=cell(length(s),1);

deltas=(sr2-2*s)/2;
for i=1:length(s)
   dims{i}=2+deltas(i):sr2(i)-deltas(i);
end

gh11=gh11(dims{:});
nh11=nh11(dims{:});

end