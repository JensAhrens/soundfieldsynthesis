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

% Fig. 5.9(a,b) and 5.10(a,b) will be created; refer also to the comment on
% lines 30-31

f     = 1000;
omega = 2*pi*f;
c     = 343;
y_ref = 1;
y_s   = -1;
x_s   =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nModify line 32 if you run out of memory.\n' );

% Make the spatial interval larger for a smoother result; make it smaller
% for lower memory requirements
spatial_interval = [ -300 300 ];

delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval(1) : delta_x : spatial_interval(2);
Y = linspace( -1, 3, 201 );

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

% initialize D_kx
D_kx = zeros( size( k_x ) );

% Eq. (4.52), (4.53), first case
D_kx( abs( k_x ) <= omega/c ) = exp( 1i .* k_x( abs( k_x ) <= omega/c ) .* x_s ) .* ...
    besselh( 0, 2, sqrt( (omega/c).^2 - k_x( abs( k_x ) <= omega/c ).^2 ) * ( y_ref - y_s ) ) ./ ...
        besselh( 0, 2, sqrt( ( omega/c ).^2 - k_x( abs( k_x ) <= omega/c ).^2 ) * y_ref );

% Eq. (4.52), (4.53), second case
D_kx( abs( k_x ) > omega/c ) = exp( 1i .* k_x( abs( k_x ) > omega/c ) .* x_s ) .* ...
    besselk( 0, sqrt( k_x( abs( k_x ) > omega/c ).^2 - ( omega/c ).^2 ) * ( y_ref - y_s ) ) ./ ...
       besselk( 0, sqrt( k_x( abs( k_x ) > omega/c ).^2 - ( omega/c ).^2 ) * y_ref );
    
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

% Fig. 5.9(a)
figure;
imagesc( X, Y, real( S ), [ -2 2 ] )

% plot secondary source distribution
hold on;
plot( [ -2 2 ], [ 0 0 ], 'k', 'LineWidth', 2 );
hold off;

colormap gray;
xlim( [ -2 2 ] );
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

graph_defaults;

% Fig. 5.9(b)
figure;
imagesc( X, Y, 20*log10( abs( S ) ), [ -10 10 ] );
turn_imagesc;
colormap gray;
revert_colormap;
xlim( [ -2 2 ] );
colorbar;
axis square;

% plot secondary source distribution
hold on;
plot( [ -2 2 ], [ 0 0 ], 'k', 'LineWidth', 2 );
hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

% Fig. 5.10(a)
figure;
plot( X, 20*log10( abs( S( (end+1)/2, : ) ) ), 'LineWidth', 2, 'Color', [ .7 .7 .7 ] ); 
turn_imagesc;
colormap gray;
revert_colormap;
axis square;

% plot desired sound field
hold on;
r     = sqrt( ( X - x_s ).^2 + ( y_ref - y_s ).^2 );
S_des = exp( -1i .* omega/c .* r ) ./ ( r + eps ); % avoid division by 0
S_des = S_des ./ S_des( (end+1)/2 ); % normalization
plot( X, 20*log10( abs( S_des ) ), 'k--' );
hold off;

grid on;

ylim( [ -6 4 ] );
xlim( [ -2 2 ] );

xlabel( 'x (m)' );
graph_defaults;

% Fig. 5.10(b)
figure;
plot( Y, 20*log10( abs( S( :, (end+1)/2 ) ) ), 'LineWidth', 2, 'Color', [ .7 .7 .7 ] ); 
turn_imagesc;
colormap gray;
revert_colormap;
axis square;

% plot desired sound field
hold on;
r     = sqrt( ( 0 - x_s ).^2 + ( Y - y_s ).^2 );
S_des = exp( -i .* omega/c .* r ) ./ ( r + eps ); % avoid division by 0
S_des = S_des ./ S_des( (end+1)/2 ); % normalization
plot( Y, 20*log10( abs( S_des ) ), 'k--' );
hold off;

grid on;

ylim( [ -10 20 ] );
xlim( [ -1   3 ] );
xlabel( 'y (m)' );
graph_defaults;
