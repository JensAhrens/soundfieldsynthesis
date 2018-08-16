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

% Note that you can create your own complex source by choosing different 
% eta and weight. They can be of arbitrary length, which has to equal for 
% both quantities though.

eta    = [ 1 4 ]; 
weight = [ 1 .5 + 1i*1.5 ];

f     = 1000;
omega = 2*pi*f;
c     = 343;
k_pw  = omega/c;
L     = 4;
y_s   = -1;
y_ref =  1;

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -10 10 ];
delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
Y = linspace( -1, 5, 401 );

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

%%%%%%%%%%%%%%%%%%% calculate desired sound field %%%%%%%%%%%%%%%%%%%%%%%%%
% initialize w_tilde
w_tilde = zeros( size( k_x ) );

% loop over all vibration modes
for index = 1 : length( eta )

    % loop over vibrating sections
    for l = 0 : eta( index )

        % Eq. (5.21)
        x_2 = -L/2 + (l+1)*L / ( eta( index ) + 1 );
        x_1 = -L/2 + (l  )*L / ( eta( index ) + 1 );
                                
        % Eq. (5.24) with first case of (5.25)
        w_tilde( k_x ~= 0 ) = w_tilde( k_x ~= 0 ) + weight( index ) .* (-1)^l .* ...
                  ( exp( 1i .* k_x( k_x ~= 0 ) .* x_2 ) - exp( 1i .* k_x( k_x ~= 0 ) .* x_1 ) ) ./ ...
                       ( 1i .* k_x( k_x ~= 0 ) );

        % Eq. (5.24) with second case of (5.25)
        w_tilde( k_x == 0 ) = w_tilde( k_x == 0 ) + weight( index ) .* (-1)^l .* ( x_2 - x_1 );
    end

end

% initialize G_kx
G_kx = zeros( size( k_x ) );

% Eq. (C.10), first case
G_kx( abs( k_x ) <= omega/c ) = -1i/4 * ...
    besselh( 0, 2, sqrt( ( omega/c ).^2 - k_x( abs( k_x ) <= omega/c ).^2 ) .* ( y_ref + abs( y_s ) ) );

% Eq. (C.10), second case
G_kx( abs( k_x ) > omega/c ) = 1/(2*pi) * ...
    besselk( 0, sqrt( k_x( abs( k_x ) > omega/c ).^2 - ( omega/c ).^2 ) .* ( y_ref + abs( y_s ) ) );

% Eq. (5.24) and (5.26), respectively
S_kx = w_tilde .* G_kx; % S_kx denotes desired sound field along y_ref
%%%%%%%%%%%%%%%%%%% end calculate desired sound field %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% calculate synthesized sound field %%%%%%%%%%%%%%%%%%%%%
% initialize D_kx
D_kx = zeros( size( k_x ) );

% Eq. (3.72) (for omnidirectional secondary source, Eq. (C.10), first case)
D_kx( abs( k_x ) <= omega/c ) = S_kx( abs( k_x ) <= omega/c ) ./ ...
    ( -1i/4 * besselh( 0, 2, sqrt( ( omega/c ).^2 - k_x( abs( k_x ) <= omega/c ).^2 ) * y_ref ) );
  
% Eq. (3.72) (for omnidirectional secondary source, Eq. (C.10), second case)
D_kx( abs( k_x ) > omega/c ) = S_kx( abs( k_x ) > omega/c ) ./ ...
    ( 1/(2*pi) * besselk( 0, sqrt( k_x( abs( k_x ) > omega/c ).^2 - ( omega/c ).^2 ) * y_ref ) );

% initialize G_kx
G_kx = zeros( size( y_m ) );

% Eq. (C.10), first case
G_kx( abs( k_x_m ) <= omega/c ) = -1i/4 * ...
    besselh( 0, 2, sqrt( ( omega/c ).^2 - k_x_m( abs( k_x_m ) <= omega/c ).^2 ) .* y_m( abs( k_x_m ) <= omega/c ) );

% Eq. (C.10), second case
G_kx( abs( k_x_m ) > omega/c ) = 1/(2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs( k_x_m ) > omega/c ).^2 - ( omega/c ).^2 ) .* y_m( abs( k_x_m ) > omega/c ) );

% synthesized sound field, Eq. (3.71)
S_kx = repmat( D_kx, [ size( G_kx, 1 ) 1 ] ) .* G_kx;
S    = ifftx( S_kx, [], 2 );
%%%%%%%%%%%%%%%%% end calculate synthesized sound field %%%%%%%%%%%%%%%%%%%

% normalization 
S = S ./ abs( S( (end+1)/2, (end+1)/2 ) );

% apply phase shift of 45 deg to make it look exactly like Fig. 5.26
S = S .* sqrt( -1i );

figure;
imagesc( X, Y, real( S ), [ -1 1 ] )
colormap gray;
xlim( [ -3 3 ] );
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

hold on;
% plot secondary source distribution
plot( [ -3 3 ], [ 0 0 ], 'k', 'LineWidth', 2 );
hold off;

graph_defaults;
