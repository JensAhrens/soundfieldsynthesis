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

clear;

f     =  700;% for Fig. 4.41(a)
%f     = 3000;% for Fig. 4.41(b)

omega = 2*pi*f;
c     = 343;
y_ref = 1;
y_s   = -1;
x_s   =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -50 50 ];
delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
Y = linspace( -2, 2, 201 );

k_x_s = (2*pi) / delta_x; % spatial sampling frequency

% create k_x
k_x    = linspace( 0, k_x_s/2, ( length( X ) + 1 ) / 2 ); % positive frequencies
k_x(1) = k_x(2); % to avoid numerical instabilities
k_x    = [ -fliplr( k_x( 2 : end ) ), k_x ]; % adds negative frequencies

% create 2D grid
[ k_x_m, y_m ]    = meshgrid( k_x, Y );
y_m               = abs( y_m );
y_m( y_m < 0.01 ) = 0.01; % to avoid numerical instablilities
%%%%%%%%%%%%%%%%%%%%%%% end prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize G_kx
G_kx = zeros( size( y_m ) );

% Eq. (C.10), first case
G_kx( abs( k_x_m ) <= omega/c ) = -1i/4 * ...
    besselh( 0, 2, sqrt( (omega/c).^2 - k_x_m( abs( k_x_m ) <= omega/c ).^2 ) .* y_m( abs( k_x_m ) <= omega/c ) );

% Eq. (C.10), second case
G_kx( abs( k_x_m ) > omega/c ) = 1/(2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs( k_x_m ) > omega/c ).^2 - (omega/c).^2 ) .* y_m( abs( k_x_m ) > omega/c ) );

% Eq. (4.55)
G_kx( abs(k_x_m) > 15 ) = 0;

G = ifftx( G_kx, [], 2 );

% normalization
G = G ./ abs( G( (end-1)/4, (end+1)/2 ) );

figure;
imagesc( X, Y, real( G ), [ -2 2 ] );

hold on;
% plot position of source
plot( 0, 0, 'kx' );
hold off;

colormap gray;
xlim( [ -2 2 ] );
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

graph_defaults;
