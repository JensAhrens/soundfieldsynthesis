%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

subfigure = 'a'; % for Fig. 5.29(a)
%subfigure = 'b'; % for Fig. 5.29(b)

f     = 1000;
omega = 2*pi*f;
c     = 343;
d_ref = 1;
y_s   = -1; % 1 m behind the secondary source distribution
x_s   =  0;

N = 21;
alpha_o = pi/3;
beta_o  = pi/2;
alpha_n = pi/2;

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -100 100 ];

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
alpha_0 = atan2( y_0, X );

% initialize D
D = zeros( size( X ) );

% this loop evaluates Eq. (5.18)
for n = 0 : N-1
    
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    for m = -n : n

        % Eq. (2.44)
        S_breve = 1i^( -n ) .* ...
                ( factorial( N - 1 ) * factorial( N ) ) ./ ...
                    ( factorial( N + n ) * factorial( N - n - 1 ) ) .* ...
                        sphharm( n, -m, beta_o, alpha_o );

        % Eq. (5.18)
        D = D + sqrt( (2*pi*d_ref) / ( 1i*omega ) ) .* ...
                    cos( alpha_n - alpha_0 ) .* exp( -1i .* omega/c .* r_0 ) ./ r_0 .* ...
                        1i^n .* S_breve .* sphharm( n, m, pi/2, alpha_0 );
    end
    
end

if ( subfigure == 'b' )
    % Eq. (5.37)
    D = conj( D );
end
    
% We go via Eq. (3.71) to be able to handle a continuous secondary source
% distribution.
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
S    = S ./ S( (end+1)/2, (end+1)/2 );

if ( subfigure == 'a' )
    S = S ./ 4;
elseif ( subfigure == 'b' )
    S = S .* 2;
end 

figure;
imagesc( X, Y - y_s, real( S ), [ -2 2 ] )

hold on;

% plot secondary source distribution
plot( [ -2 2 ], [ 1 1 ], 'k', 'LineWidth', 2 );

if ( subfigure == 'b' )
    % plot boundary between converging and diverging fields
    plot( [ -2 2 ], [ 2 2 ], 'k:' );
    
    % plot focus point
    plot( 0, 2, 'kx' );
end

hold off;

colormap gray;
xlim( [ -2 2 ] );
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

graph_defaults;
