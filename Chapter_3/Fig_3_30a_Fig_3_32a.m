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

% Note that Fig. 3.30(a) and 3.32(a) are identical.

clear;

N = 60; % order
f = 1200;
A = .3;
c = 343;
k = 2*pi*f/c;

theta_pw = pi/2;
phi_pw   = pi/2;

% create spatial grid
resolution = 200;
X          = linspace( -2, 2, resolution );
Y          = linspace( -2, 2, resolution );
[ x, y ]   = meshgrid( X, Y );
r          = sqrt( x.^2 + y.^2 );
alpha      = atan2( y, x );
beta       = pi/2;

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
        S_breve = 4 .* pi .* 1i^(-n) .* sphharm( n, -m, phi_pw, theta_pw );

        % Eq. (3.108)
        S_breve_scattered = - bessel_prime ./ hankel_prime .* S_breve;
        
        % this combines Eq. (2.32a) and (2.32b) 
        % ( incoming sound field plus outgoing scattered sound field )
        S = S + ( S_breve .* bessel + S_breve_scattered .* hankel ) .* sphharm( n, m, beta, alpha ); 
    end
end

% no sound field inside the sphere
S( r < A ) = 0;

figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
colormap gray;
turn_imagesc;
axis square;

hold on;

% plot sphere
x_circ = linspace( -A, A, 100 );

plot( x_circ,  sqrt( A.^2 - x_circ.^2 ), 'Color', [ .99 .99 .99 ], 'LineWidth', 2 );
plot( x_circ, -sqrt( A.^2 - x_circ.^2 ), 'Color', [ .99 .99 .99 ], 'LineWidth', 2 );

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;
