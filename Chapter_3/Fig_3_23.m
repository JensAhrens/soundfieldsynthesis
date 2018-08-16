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

% Note that this script is actually a simplified version of what is shown
% in Fig. 3.23. It is the 2D equivalent published in
%
% Jens Ahrens and Sascha Spors: On the Secondary Source Type Mismatch in
% Wave Field Synthesis Employing Circular Distributions of Loudspeakers.
% 127th Convention of the AES, New York, NY, USA, Oct. 2009

subfigure = 'a'; % for Fig. 3.23(a)
%subfigure = 'b'; % for Fig. 3.23(b)

theta_pw = pi/2;
M        = 40;
R        = 1.5;  % m
c        = 343;  % m/s

if ( subfigure == 'a' )
    f    = 200; % Hz
elseif ( subfigure == 'b' )
    f    = 1000; % Hz
end

k = (2*pi*f)/c;

alpha_0 = linspace( 0, 2*pi, 301 );
alpha_0 = alpha_0( 1 : end - 1 );

D_wfs = -2 .* 1i .* k .* cos( theta_pw - alpha_0 ) .* ...
                                exp( -1i .* k .* R .* cos( theta_pw - alpha_0 ) );
% secondary source selection
D_wfs( 1 : end/2 + 1 ) = 0;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x );

% initialize S_int and S_ext
S_int = zeros( size( r ) );
S_ext = zeros( size( r ) );

for m = -M : M
    disp( [ 'Calculating order ' num2str( m ) ' of ' num2str( M ) '.' ] );
       
    % numerical integration
    D_ring = sum( D_wfs .* exp( -1i .* m .* alpha_0 ), 2 ) / ( length( alpha_0 ) / 2 );
    
    G_ring_int = besselh( m, 2, k*R );
    G_ring_ext = besselj( m,    k*R );

    S_int = S_int + D_ring .* G_ring_int .* besselj( m,    k.*r ) .* exp( 1i * m * alpha );
    S_ext = S_ext + D_ring .* G_ring_ext .* besselh( m, 2, k.*r ) .* exp( 1i * m * alpha );  
                    
end

% assemble exterior and interior field
S          = zeros( size( x ) ); 
S( r < R ) = S_int( r < R );
S( r > R ) = S_ext( r > R );


% rotate the phase a bit to emphasize the deviation in Fig. 3.23(a).
S = S .* exp( 1i * 2*pi *.8 );

figure;
imagesc( X, Y, real( S ), [ -1 1 ] );
turn_imagesc;
colormap gray;
axis square;

hold on;

% plot secondary source distribution
x_circ = linspace( -R, R, 200 );
plot( x_circ,  sqrt( R.^2 - x_circ.^2 ), 'k:', 'Linewidth', 2 )
plot( x_circ, -sqrt( R.^2 - x_circ.^2 ), 'k' , 'Linewidth', 2 )

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;
