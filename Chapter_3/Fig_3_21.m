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

subfigure = 'a'; % for Fig. 3.21(a)
%subfigure = 'b'; % for Fig. 3.21(b)

fs    = 22050; % sampling time-frequency
c     = 343;
y_ref = 1;

L = 2;
theta_pw = pi/4;
phi_pw   = pi/2;

xlims = [ -1.5, -.5 ]; % for plotting

no_of_bins = 101 ; % half number of time-frequency sampling points

f      = linspace( 0, fs/2, no_of_bins );
f( 1 ) = f( 2 ); % to avoid numerical instabilities

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -20 20 ];
delta_x          = .005; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
Y = linspace( -.2, .8, 201 );

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

if ( subfigure == 'b' )
    w = [ zeros( 1, length( spatial_interval(1) : delta_x : -L/2 ) ) ...
          sin( linspace( 0, pi/2, length( -L/2 + delta_x : delta_x : -3*L/8 ) ) ).^2 ...
          ones( 1, length( -3*L/8 + delta_x : delta_x : 3*L/8-delta_x ) ) ...
          cos( linspace( 0, pi/2, length( 3*L/8 : delta_x : L/2 - delta_x ) ) ).^2 ...
          zeros( 1, length( L/2 : delta_x : spatial_interval(2) ) ) ];

    w_tilde = fftx( w, [] , 2 );

end

% initialize
S = zeros( length( Y ), ( xlims(2) - xlims(1) )./delta_x + 1, length( f ) ); % Y, X, f

% loop over frequencies
for index = 2 : length( f )
    disp( [ 'Calculating frequency bin ' num2str( index ) ' of ' num2str( length( f ) ) '.' ] );
    
    omega = 2*pi*f(index);
    
    k_pw_x = omega./c .* cos( theta_pw ) .* sin( phi_pw );
    k_pw_y = omega./c .* sin( theta_pw ) .* sin( phi_pw );

    D_kx = zeros( size( k_x ) );

    % no tapering
    if ( subfigure == 'a' )
        % this mimicks convolution of Eq. (3.74) with (3.87); as it is
        % stated in (3.85)
        D_kx = L .* sinc( ( k_x - k_pw_x ) .* L ./ (2*pi) ) .* ...
                        4 .* 1i .* exp( -1i .* k_pw_y .* y_ref ) ./ ...
                                            besselh( 0, 2, k_pw_y.*y_ref );
    
    elseif ( subfigure == 'b' )
        % this mimicks convolution of w_tilde with (3.87); as it is
        % stated in (3.85); the circsift represents the shift due to
        % convolution with a Dirac
        D_kx = circshift( w_tilde, [ 0 round( k_pw_x / d_k_x ) ] ) .* ...
                    4 .* 1i .* exp( -1i .* k_pw_y .* y_ref ) ./ besselh( 0, 2, k_pw_y .* y_ref );
    end

    G_kx = zeros( size( y_m ) );

    % Eq. (C.10), first case
    G_kx( abs( k_x_m ) < omega/c ) = -1i/4 * ...
        besselh( 0, 2, sqrt( ( omega/c ).^2 - k_x_m( abs( k_x_m ) < omega/c ).^2 ) .* y_m( abs( k_x_m ) < omega/c ) );
    
    % Eq. (C.10), second case
    G_kx( abs( k_x_m ) > omega/c ) = 1 / (2*pi) * ...
        besselk( 0, sqrt( k_x_m( abs( k_x_m ) > omega/c ).^2 - ( omega/c ).^2 ) .* y_m( abs( k_x_m ) > omega/c ) );

    % Eq. (3.71)
    S_kx = repmat( D_kx, [ size( G_kx, 1 ) 1 ] ) .* G_kx;
    
    clear D_kx G_kx; 
    
    S_tmp = ifftx( S_kx, [], 2 );

    % save relevant portion of result
    S( :, :, index ) = S_tmp( :, ( end+1 ) / 2 + xlims(1) / delta_x : ( end+1 ) / 2 + xlims(2) / delta_x );
    
    clear S_tmp; 
    
end

clear k_x k_x_m y_m X;

% create time-domain signal; 22 samples equals 1 ms, which equals 30 cm
chunk_length  = 2 * no_of_bins - 2;
% concatenate 10 chunks
signal        = repmat( [ sin( linspace( 0, pi, 5 ) ).^2 zeros( 1, 15 ) ].', [ chunk_length/20 1 ] );

signal_fft  = fft( signal );
signal_spec = signal_fft( 1 : end/2+1 );

% perform fast convolution of signal and sound field; watch out for 
% time-domain aliasing 
disp( 'Performing convolution.' );

for n = 1 : size( S, 1 )
    for m = 1 : size( S, 2 )
        S( n, m, : ) = squeeze( S( n, m, : ) ) .* signal_spec;
    end
end

% transform sound field to time domain
s = ifft( cat( 3, S, conj( flipdim( S( :, :, 2 : end-1 ), 3 ) ) ), [], 3 );

% normalize
if ( subfigure == 'a' )
    s = s ./ 2.5;
elseif ( subfigure == 'b' )
    s = s .* 2.5;
end

% time instant to plot in samples
time = 8;

figure;
X_plot = xlims( 1 ) : delta_x : xlims( 2 );
imagesc( X_plot, Y, real( s( :, :, time ) ), [ -2 2 ] ) 
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
