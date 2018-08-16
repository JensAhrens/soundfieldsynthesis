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

subfigure = 'a'; % for Fig. 4.36(a)
%subfigure = 'b'; % for Fig. 4.36(b)
%subfigure = 'c'; % for Fig. 4.36(c)
%subfigure = 'd'; % for Fig. 4.36(d)

f     = 1300;
omega = 2*pi*f;
c     = 343;
y_ref = 1;
y_s   = -1;
x_s   =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -75 75 ];
delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
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
    

%%%%%%%%%%%%%%%%%% bandwidth limitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( subfigure == 'c' )
    B_low  = -180;
    B_high =  180;
elseif ( subfigure == 'd' )    
    B_low  =  0;
    B_high = 300;
end
    
if ( subfigure == 'c' || subfigure == 'd' )    
    central_bin = round( length( D_kx ) / 2 );

    indices = [ 1 : central_bin + B_low , central_bin + B_high : length( D_kx ) ];
    D_kx( :, indices ) = 0;

end
%%%%%%%%%%%%%%%%%% end bandwidth limitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( subfigure == 'b' || subfigure == 'c' || subfigure == 'd' )    

    D = ifftx( D_kx, [], 2 );

    indices = 1 : length( D );
    indices = indices( mod( indices-1, 20 )~=0 );
    
    % use every 20th value; this correponds to a sampling step of 20 cm 
    D( :, indices ) = 0;
    
    D_kx = fftx( D, [], 2 );

end
%%%%%%%%%%%%%%%%%%%%%%%% end discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
S    = S ./ abs( S( (end+1)/2, (end+1)/2 ) ) .* 1.3;

figure;
imagesc( X, Y, real( S ), [ -3 3 ] )

hold on;

% plot secondary source distribution
if ( subfigure == 'a')
    plot( [ -2 2 ], [ 0 0 ], 'k', 'LineWidth', 2 );
else
    plot( -2 : 0.2 : 2, 0, 'kx', 'Color', [ .1 .1 .1 ] );
end

hold off;

colormap gray;
xlim( [ -2 2 ] );
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

graph_defaults;
