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
% (c) 2012 by Jens Ahrens                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

theta_pw = pi/2;
phi_pw   = pi/2;
N        = 50;
R        = 1.5;  % m
A        = .3;
f        = 1200; % Hz
c        = 343;  % m/s

k = (2*pi*f)/c;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x );
beta  = pi/2;

% initialize S
S = zeros( size( r ) );

for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
   
    bessel = sphbesselj( n,    k.*r );
    hankel = sphbesselh( n, 2, k.*r );
    
    % Eq. (2.18)
    bessel_prime = 1/(2*n+1) * ( n * sphbesselj( n-1,    k*A ) - (n+1) * sphbesselj( n+1,    k*A ) ) * k;
    hankel_prime = 1/(2*n+1) * ( n * sphbesselh( n-1, 2, k*A ) - (n+1) * sphbesselh( n+1, 2, k*A ) ) * k;
    
    for m = -n : n
        
        % Eq. (2.38)
        S_breve = 4 .* pi .* 1i^(-n) .* sphharm( n, -m, pi/2, theta_pw );

        % synthesized interior field, no scatterer, Eq. (3.22) 
        S( r<R ) = S( r<R ) + S_breve .* bessel( r<R ) .* sphharm( n, m, beta, alpha( r<R ) ); 
       
        % from Eq. (2.37a)
        G_breve_int = -1i .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, 0, 0, 0 );
       
        % from Eq. (2.37b)
        G_breve_ext = -1i .* k .* sphbesselj( n, k.*R )    .* sphharm( n, 0, 0, 0 );
       
        % synthesized exterior field, no scatterer, Eq. (3.23)
        S( r>R ) = S( r>R ) + S_breve .* G_breve_ext ./ G_breve_int .* hankel( r>R ) .* sphharm( n, m, beta, alpha( r>R ) );
        
        % Eq. (3.112)
        S_breve_scattered = - bessel_prime ./ hankel_prime .* S_breve;
        
        % scattered field
        S( r>A ) = S( r>A ) + S_breve_scattered .* hankel( r>A ) .* sphharm( n, m, beta, alpha( r>A ) );
        
    end
    
end

% no sound field inside the sphere
S( r < A ) = 0;

figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
turn_imagesc;
colormap gray;
axis square;

hold on;

% plot secondary source distribution
x_circ = linspace( -R, R, 200 );
plot( x_circ,  sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )
plot( x_circ, -sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )

% plot sphere
x_circ = linspace( -A, A, 100 );

plot( x_circ,  sqrt( A.^2 - x_circ.^2 ), 'Color', [ .99 .99 .99 ], 'LineWidth', 2 );
plot( x_circ, -sqrt( A.^2 - x_circ.^2 ), 'Color', [ .99 .99 .99 ], 'LineWidth', 2 );

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;
