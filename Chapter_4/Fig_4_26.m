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

subfigure = 'a'; % for Fig. 4.26(a)
%subfigure = 'b'; % for Fig. 4.26(b)

if ( subfigure == 'a' )
    f = 1000; 
else
    f = 3000;
end

theta_pw = pi/2;
phi_pw   = pi/2;
N        = 79;
M        = 27;
R        = 1.5;  % m
L        = 56;
c        = 343;  % m/s

k = (2*pi*f)/c;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2);
alpha = atan2( y, x );
beta  = pi/2;

% initialize S_int and S_ext
S_int = zeros( size( r ) );
S_ext = zeros( size( r ) );

% this loop evaluates Eq. (3.22) and (3.23)
for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
   
    for m = max( -n, -M ) : min( n, M )
       
        % Eq. (2.38)
        S_breve = 4*pi .* 1i^( -n ) .* sphharm( n, -m, phi_pw, theta_pw );
       
        S_int = S_int + S_breve .* sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );

        % from Eq. (2.37a)
        G_breve_int = -i .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, 0, 0, 0 );
       
        % from Eq. (2.37b)
        G_breve_ext = -1i .* k .* sphbesselj( n, k.*R )    .* sphharm( n, 0, 0, 0 );
       
        S_ext = S_ext + S_breve .* G_breve_ext ./ G_breve_int .* ...
                             sphbesselh( n, 2, k.*r ) .* sphharm( n, m, beta, alpha );
        
    end
end

% assemble exterior and interior field
S          = zeros( size( x ) ); 
S( r < R ) = S_int( r < R );
S( r > R ) = S_ext( r > R );

figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
turn_imagesc;
colormap gray;
axis square;

hold on;

% plot secondary source distribution
alpha_0 = linspace( 0, 2*pi, L+1 );
alpha_0 = alpha_0( 1 : end-1 );

plot( R*cos( alpha_0 ),  R*sin( alpha_0 ), 'kx' )

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;
