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

% both Fig. 3.19(a) and 3.19(b) will be created

f     = 1000;
c     = 343;
y_ref = 1;
L     = 2;

theta_pw = pi/4;
phi_pw   = pi/2;

k = 2*pi*f/c;
    
%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -120 120 ];
delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
Y = linspace( -1, 5, 201 );

k_x_s = (2*pi) / delta_x; % spatial sampling frequency

d_k_x = k_x_s / length( X ); % sampling interval

% create k_x
k_x    = linspace( 0, k_x_s/2, ( length( X ) + 1 ) / 2 ); % positive frequencies
k_x(1) = k_x(2); % to avoid numerical instabilities
k_x    = [ -fliplr( k_x( 2 : end ) ), k_x ]; % adds negative frequencies

% create 2D grid
[ k_x_m, y_m ]    = meshgrid( k_x, Y );
y_m               = abs( y_m );
y_m( y_m < 0.01 ) = 0.01; % to avoid numerical instablilities
%%%%%%%%%%%%%%%%%%%%%%% end prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize S
S = zeros( length( Y ), length( X ) );

k_pw_x = k .* cos( theta_pw ) .* sin( phi_pw );
k_pw_y = k .* sin( theta_pw ) .* sin( phi_pw );

D_kx = zeros( size( k_x ) );

% this mimicks convolution of Eq. (3.74) with (3.87); as it is
% stated in (3.85)
D_kx = L .* sinc( ( k_x - k_pw_x ) .* L ./ (2*pi) ) .* ...
                4.*1i .* exp( -1i .* k_pw_y .* y_ref ) ./ ...
                                    besselh( 0, 2, k_pw_y.*y_ref );

% omnidirectional secondary source
G_kx = zeros( size( y_m ) );

% Eq. (C.10), first case
G_kx( abs( k_x_m ) < k ) = -1i/4 * ...
    besselh( 0, 2, sqrt( k.^2 - k_x_m( abs(k_x_m) < k ).^2 ) .* y_m( abs(k_x_m) < k ) );

% Eq. (C.10), second case
G_kx( abs(k_x_m) > k ) = 1/(2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs(k_x_m) > k ).^2 - k.^2 ) .* y_m( abs(k_x_m) > k ) );

% Eq. (3.71)
S_kx = repmat( D_kx, [ size( G_kx, 1 ) 1 ] ) .* G_kx;


S = ifftx( S_kx, [], 2);

% normalize
S = S ./ abs( S( 68, find( X == 1 ) ) ); % corresponds to location x = 1, y = 1

% Fig. 3.19(a)
figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
xlim( [ -2 4 ] )
colormap gray;
revert_colormap;
colorbar;
    
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

hold on;
% plot secondary source distribution
plot( [ -L/2 L/2 ], [ 0 0 ], 'k', 'Linewidth', 2 );
hold off;

graph_defaults;

% Fig. 3.19(b)
figure;
imagesc( X, Y, 20*log10( abs( S ) ), [ -40 10 ] );
xlim( [ -2 4 ] )
colormap gray;
revert_colormap;
colorbar;

xlabel('x (m)');
ylabel('y (m)');
turn_imagesc;
axis square;

hold on;
% plot secondary source distribution
plot( [ -L/2 L/2 ], [ 0 0 ], 'k', 'Linewidth', 2 );
hold off;

graph_defaults;
