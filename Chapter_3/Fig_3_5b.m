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

theta_pw = pi/2;
phi_pw   = pi/2;
N        = 40;
N_prime  = 13;
R        = 1.5;  % m
f        = 1000; % Hz
c        = 343;  % m/s

% orientation of secondary source at pole
alpha_or = 0;
beta_or  = pi;

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

% get coefficients (E|I)_n,n'
EI = zonal_translation_coefficients( N, N_prime, k*R, 'EI' );

% get coefficients (E|E)_n,n'
EE = zonal_translation_coefficients( N, N_prime, k*R, 'EE' );

% calculate directivity coefficients of secondary sources
G_breve_ext_prime = zeros( N_prime, 1 );

% this loop evaluates Eq. (2.44) for m = 0
for n_prime = 0 : N_prime - 1     
        G_breve_ext_prime( n_prime+1 ) = ...
            1i^( -n_prime ) .* ...
                ( factorial( N_prime - 1 ) * factorial( N_prime ) ) ./ ...
                ( factorial( N_prime + n_prime ) * factorial( N_prime - n_prime - 1 ) ) .* ...
                                   sphharm( n_prime, 0, beta_or, alpha_or );
               
end

% this loop evaluates Eq. (3.22) and (3.23)
for n = 0 : N - 1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N - 1 ) '.' ] );
   
    for m = -n : n
       
        % Eq. (2.38)
        S_breve = 4*pi .* 1i^( -n ) .* sphharm( n, -m, phi_pw, theta_pw );
        
        % Eq. (3.22) 
        S_int = S_int + S_breve .* sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );
                
        % initialize G_breve_int and G_breve_ext
        G_breve_int = 0;
        G_breve_ext = 0;
        
        % this loop evaluates Eq. (3.27), see also the erratum at
        % http://soundfieldsynthesis.org/?page_id=84
        for n_prime = 0 : N_prime - 1
            
            G_breve_int = G_breve_int + G_breve_ext_prime( n_prime + 1 ) * (-1)^( n + n_prime ) * ... 
                                        EI( n + 1, n_prime + 1 );
       
            G_breve_ext = G_breve_ext + G_breve_ext_prime( n_prime + 1 ) * (-1)^( n + n_prime ) * ... 
                                        EE( n + 1, n_prime + 1 );
       
        end
        
        % Eq. (3.23)
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
x_circ = linspace( -R, R, 200 );
plot( x_circ,  sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )
plot( x_circ, -sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );

graph_defaults;
