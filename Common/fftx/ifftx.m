function [ out ] = ifftx( varargin )
%IFFTX Spatial fourier transform with negative exponent; 
%    it assumes that the zero of the space-frequency axis is in the
%    center of the buffer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work is supplementary material for the book                        %
%                                                                         %
% Jens Ahrens, Analytic Methods of Sound Field Synthesis, Springer-Verlag %
% Berlin Heidelberg, 2012, http://dx.doi.org/10.1007/978-3-642-25743-8    %
%                                                                         %
% It has been downloaded from http://soundfieldsynthesis.org and is       %
% licensed under a Creative Commons Attribution-NonCommercial-ShareAlike  % 
% 3.0 Unported License. Please cite the book appropriately if you use     % 
% these materials in your own work.                                       %
%                                                                         %
% (c) 2012 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    input_data = varargin{1};
    points     = [];
    dimension  = 1;

elseif nargin == 2
    input_data = varargin{1};
    points     = varargin{2};
    dimension  = 1;
    
elseif nargin == 3
    input_data = varargin{1};
    points     = varargin{2};
    dimension  = varargin{3};
    
else
    error( 'Too many input arguments.' );
end
    
% shift zero-frequency to first bin as fft expects it 
shift            = zeros( 1, ndims(input_data) );
shift(dimension) = ceil( size( input_data, dimension ) / 2 );
input_data       = circshift( input_data, shift );

% do FFT
out = fft( input_data, points, dimension );

% put zero of space axis to the center of the buffer
out = circshift( out, -shift );

% normalization
out = out ./ (2*pi);

end
