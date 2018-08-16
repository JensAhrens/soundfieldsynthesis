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

% Note that this script uses the myspectrogram function from 
% https://ccrma.stanford.edu/~jos/sasp/Matlab_listing_myspectrogram_m.html

clear;

v     = 40;
f     = 500; 
fs    = 44100;
c     = 343; 
omega = 2*pi*f;
M     = v / c;

x = 0;
y = 4;
z = 0;

t = -1 : 1 / fs : 1;

% Eq. (5.58)
Delta = sqrt( ( x - v*t ).^2  + ( y.^2 + z.^2 ) .* ( 1 - M^2 ) );

% Eq. (5.57)
tau = ( M .* ( x - v*t ) + Delta ) / ( c * ( 1 - M^2 ) );

% Eq. (5.60)
s = 1 / (4*pi) .* exp( 1i .* omega .* ( t - tau ) ) ./ Delta;

% normalize
s = s ./ max( abs( s ) );

win_length = round( fs / 5 ); 
hop        = round( win_length / 32 ); 

spec      = myspectrogram( real( s ), win_length, fs, hamming( win_length ), win_length - hop, false );
size_spec = size( spec );

t_plot = linspace( t( 1 ), t( end ), size_spec( 2 ) );
f_plot = linspace(      0,       fs, size_spec( 1 ) );

figure;

imagesc( t_plot, f_plot, 20*log10( abs( spec ) / win_length ), [ -60 -20 ] );

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
grid on;
axis square;

xlim(  [ -0.5 0.5 ] );
ylim(  [  350 650 ] );

xlabel( 't (s)' );
ylabel( 'f (Hz)' );

graph_defaults;
