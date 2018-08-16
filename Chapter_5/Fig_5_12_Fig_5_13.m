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

% Fig. 5.12(a,b) and 5.13(a,b) will be created; refer also to the comment on
% lines 30-31

f     = 1000;
omega = 2*pi*f;
c     = 343;
d_ref = 1;
y_s   = -1; % 1 m behind the secondary source distribution
x_s   =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nModify line 32 if you run out of memory.\n' );

% Make the spatial interval larger for a smoother result; make it smaller
% for lower memory requirements
spatial_interval = [ -300 300 ];

delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );

% we calculate as if the secondary source distribution were along the x
% axis and the source at y = -1
Y = linspace( -1, 3, 201 ); 

k_x_s = (2*pi) / delta_x; % spatial sampling frequency

% create k_x
k_x    = linspace( 0, k_x_s/2, ( length( X ) + 1 ) / 2 ); % positive frequencies
k_x(1) = k_x(2); % to avoid numerical instabilities
k_x    = [ -fliplr( k_x( 2 : end ) ), k_x ]; % adds negative frequencies

% create 2D grids
[ k_x_m, y_m ]    = meshgrid( k_x, Y );
y_m               = abs( y_m );
y_m( y_m < 0.01 ) = 0.01; % to avoid numerical instablilities
%%%%%%%%%%%%%%%%%%%%%%% end prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%

y_0 = -y_s;
r_0 = sqrt( X.^2 + y_0.^2 );

% Eq. (5.14)
D = sqrt( 2 * pi * d_ref * 1i * omega ) .* y_0 ./ r_0 .* ...
                                    exp( -1i .* omega ./ c .* r_0 ) ./ r_0;  

% we go via Eq. (3.71) to be able to handle a continuous sceondary source
% distribution
D_kx = fftx( D, [], 2 );
    
% initialize G_kx
G_kx = zeros( size( y_m ) );

% Eq. (C.10), first case
G_kx( abs( k_x_m ) <= omega/c ) = -1i/4 * ...
    besselh( 0, 2, sqrt( (omega/c).^2 - k_x_m( abs( k_x_m ) <= omega/c ).^2 ) .* y_m( abs( k_x_m ) <= omega/c ) );

% Eq. (C.10), second case
G_kx( abs( k_x_m ) > omega/c ) = 1/(2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs( k_x_m ) > omega/c ).^2 - (omega/c).^2 ) .* y_m( abs( k_x_m ) > omega/c ) );

% Eq. (3.71)
S_kx = repmat( D_kx, [ size( G_kx, 1 ) 1 ] ) .* G_kx;
S    = ifftx( S_kx, [], 2 );

% normalization
S    = S ./ abs( S( (end+1)/2, (end+1)/2 ) );

% Fig. 5.12(a)
figure;
imagesc( X, Y - y_s, real( S ), [ -2 2 ] )

% plot secondary source distribution
hold on;
plot( [ -2 2 ], [ 1 1 ], 'k', 'LineWidth', 2 );
hold off;

colormap gray;
xlim( [ -2 2 ] );
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

graph_defaults;

% Fig. 5.12(b)
figure;
imagesc( X, Y - y_s, 20*log10( abs( S ) ), [ -10 10 ] );
turn_imagesc;
colormap gray;
revert_colormap;
xlim( [ -2 2 ] );
colorbar;
axis square;

% plot secondary source distribution
hold on;
plot( [ -2 2 ], [ 1 1 ], 'k', 'LineWidth', 2 );
hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

% Fig. 5.13(a)
figure;
plot( X, 20*log10( abs( S( (end+1)/2, : ) ) ), 'Linewidth', 2, 'Color', [ .7 .7 .7 ] ); 
turn_imagesc;
colormap gray;
revert_colormap;
axis square;

% plot desired sound field
hold on;
r     = sqrt( ( X - x_s ).^2 + ( d_ref - y_s ).^2 );
S_des = exp( -1i .* omega/c .* r ) ./ ( r + eps ); % avoid division by 0
S_des = S_des ./ S_des( (end+1)/2 ); % normalization
plot( X, 20*log10( abs( S_des ) ), 'k--' );
hold off;

grid on;

ylim( [ -6 4 ] );
xlim( [ -2 2 ] );

xlabel( 'x (m)' );
graph_defaults;

% Fig. 5.13(b)
figure;
plot( Y + 1, 20*log10( abs( S( :, (end+1)/2 ) ) ), 'Linewidth', 2, 'Color', [ .7 .7 .7 ] ); 
turn_imagesc;
colormap gray;
revert_colormap;
axis square;

% plot desired sound field
hold on;
r     = sqrt( ( 0 - x_s ).^2 + ( Y - y_s ).^2 );
S_des = exp( -1i .* omega/c .* r ) ./ ( r + eps ); % avoid division by 0
S_des = S_des ./ S_des( (end+1)/2 ); % normalization
plot( Y + 1, 20*log10( abs( S_des ) ), 'k--' );
hold off;

grid on;

ylim( [ -10 20 ] );
xlim( [   0  4 ] );

xlabel( 'y (m)' );
graph_defaults;
