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

subfigure = '3.26'; % for Fig. 3.26(a) and (b) 
%subfigure = '3.27'; % for Fig. 3.27(a) and (b)

f     = 1000;
c     = 343;
omega = 2*pi*f;

y_ref = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -20 20 ];
delta_x          = .005; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
Y = -.5 : delta_x : 3.5;

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

% Note that in the following it is first assumed that the secondary source
% distributions are both located along the x-axis. They are then rotated 
% appropriately.

% spatial window
if ( strcmp( subfigure, '3.26' ) )
    w = ones( 1, length( X ) );
    w( X < 0 ) = 0;
    w( X > 3 ) = 0;

elseif ( strcmp( subfigure, '3.27' ) )
    w = ones( 1, length( X ) );
    w( X < 0 )   = 0;
    w( X > 1.5 ) = cos( ( X( X > 1.5 ) - 1.5 ) / 1.5 * pi/2 ).^2;
    w( X > 3 )   = 0;
    
end

w_tilde = fftx( w, [] , 2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% portion that is parallel to x-axis %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_pw = pi/6;
phi_pw   = pi/2;

k_pw_x = omega/c * cos( theta_pw ) * sin( phi_pw );
k_pw_y = omega/c * sin( theta_pw ) * sin( phi_pw );

% this mimicks convolution of w_tilde with (3.87); as it is stated in 
% (3.85); the circshift represents the shift due to convolution with a 
% Dirac; we have enough space so that the circularity is not an issue
D_kx = circshift( w_tilde, [ 0 round( k_pw_x / d_k_x ) ] ) .* ...
            4 * 1i .* exp( -1i .* k_pw_y .* y_ref ) ./ besselh( 0, 2, k_pw_y .* y_ref );

G_kx = zeros( size( y_m ) );

% Eq. (C.10), first case
G_kx( abs( k_x_m ) < omega/c ) = - 1i / 4 * ...
    besselh( 0, 2, sqrt( ( omega/c ).^2 - k_x_m( abs(k_x_m) < omega/c ).^2 ) .* y_m( abs( k_x_m ) < omega/c ) );

% Eq. (C.10), second case
G_kx( abs( k_x_m ) > omega/c ) = 1 / (2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs( k_x_m ) > omega/c ).^2 - ( omega/c ).^2 ) .* y_m( abs( k_x_m ) > omega/c ) );

% Eq. (3.71)
S_kx = repmat( D_kx, [ size( G_kx, 1 ) 1 ] ) .* G_kx;

S_x_axis = ifftx( S_kx, [], 2 );

% remove unrequired data
S_x_axis( :, X < -.5 | X > 3.5 ) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% portion that is parallel to y-axis %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


theta_pw = 2/6 * pi;
phi_pw   = pi/2;

k_pw_x = omega/c * cos( theta_pw ) * sin( phi_pw );
k_pw_y = omega/c * sin( theta_pw ) * sin( phi_pw );

% this mimicks convolution of w_tilde with (3.87); as it is stated in 
% (3.85); the circshift represents the shift due to convolution with a 
% Dirac; we have enough space so that the circularity is not an issue
D_kx = circshift( w_tilde, [ 0 round( k_pw_x / d_k_x ) ] ) .* ...
            4 * 1i .* exp( -1i .* k_pw_y .* y_ref ) ./ besselh( 0, 2, k_pw_y .* y_ref );

G_kx = zeros( size( y_m ) );

% Eq. (C.10), first case
G_kx( abs( k_x_m ) < omega/c ) = - 1i / 4 * ...
    besselh( 0, 2, sqrt( ( omega/c ).^2 - k_x_m( abs(k_x_m) < omega/c ).^2 ) .* y_m( abs( k_x_m ) < omega/c ) );

% Eq. (C.10), second case
G_kx( abs( k_x_m ) > omega/c ) = 1 / (2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs( k_x_m ) > omega/c ).^2 - ( omega/c ).^2 ) .* y_m( abs( k_x_m ) > omega/c ) );

% Eq. (3.71)
S_kx = repmat( D_kx, [ size( G_kx, 1 ) 1 ] ) .* G_kx;

S_y_axis = ifftx( S_kx, [], 2 );

% remove unrequired data
S_y_axis( :, X < -.5 | X > 3.5 ) = [];

% rotate by 90 deg
S_y_axis = S_y_axis.';

S = S_x_axis + S_y_axis;

% subfigure (a)
figure;
imagesc( X( abs( X ) < 2 ), Y - 1.5, real( S ), [ -1 1 ] ) 
colormap gray;
revert_colormap;

xlim( [ -2 2 ] );
ylim( [ -2 2 ] );

xlabel( 'x (m)' );
ylabel( 'y (m)' );

turn_imagesc;
axis square;

hold on;
% plot secondary source distribution
plot( [ -1.5  1.5 ], [ -1.5 -1.5 ], 'k' , 'Linewidth', 2 );
plot( [ -1.5 -1.5 ], [ -1.5  1.5 ], 'k' , 'Linewidth', 2 );
plot( [ -1.5  1.5 ], [  1.5  1.5 ], 'k:', 'Linewidth', 2 );
plot( [  1.5  1.5 ], [ -1.5  1.5 ], 'k:', 'Linewidth', 2 );
hold off;

graph_defaults;

% subfigure (b)
figure;
imagesc( X( abs( X ) < 2 ), Y - 1.5, 20*log10( abs( S ) ), [ -10 10 ] ) 
colormap gray;
revert_colormap;
colorbar;

xlim( [ -2 2 ] );
ylim( [ -2 2 ] );

xlabel( 'x (m)' );
ylabel( 'y (m)' );

turn_imagesc;
axis square;

hold on;
% plot secondary source distribution
plot( [ -1.5  1.5 ], [ -1.5 -1.5 ], 'k' , 'Linewidth', 2 );
plot( [ -1.5 -1.5 ], [ -1.5  1.5 ], 'k' , 'Linewidth', 2 );
plot( [ -1.5  1.5 ], [  1.5  1.5 ], 'k:', 'Linewidth', 2 );
plot( [  1.5  1.5 ], [ -1.5  1.5 ], 'k:', 'Linewidth', 2 );
hold off;

graph_defaults;
