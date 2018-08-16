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

fprintf( ['\nWarning: This script takes awfully long to execute because of the\n' ...
          'necessity for the calculation of the general translation coefficients.\n\n' ] );

clear;

tic;

theta_pw = pi/2;
phi_pw   = pi/2;
N        = 30; 
N_prime  = 13;
R        = 1.5; % m
f        = 700; % Hz
c        = 343; % m/s

% orientation of secondary source at alpha = 0
alpha_or = pi;
beta_or  = pi/2;

k = (2*pi*f)/c;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x );
beta  = pi/2;

% get coefficients (E|I)_(n,n')^(m,m'), we need them for Eq. (3.50), (3.51)
EI = general_translation_coefficients( N, N_prime, k*R, 0, pi/2, 'EI' );

% get coefficients (E|E)_(n,n')^(m,m'), we need them for Eq. (3.50), (3.51)
EE = general_translation_coefficients( N, N_prime, k*R, 0, pi/2, 'EE' );

% get coefficients (E|E)_(|m|,n')^(m,m'), Eq. (3.55); we could also deduce 
% these coefficients from the general coefficients EI but like this it's 
% more fun...
EI_sect = sectorial_translation_coefficients( N-1, N_prime, k*R, 0, pi/2, 'EI' );

% calculate directivity coefficients of secondary sources
G_breve_prime = zeros( N_prime^2 );

% this loop evaluates Eq. (2.44)
for n_prime = 0 : N_prime - 1     
    
    for m_prime = -n_prime : n_prime

        G_breve_prime( 1 + n_prime^2 + n_prime + m_prime ) = 1i^( -n_prime ) .* ...
                ( factorial( N_prime - 1 ) * factorial( N_prime ) ) ./ ...
                ( factorial( N_prime + n_prime ) * factorial( N_prime - n_prime - 1 ) ) .* ...
                                   sphharm( n_prime, -m_prime, beta_or, alpha_or );

    end

end

toc;

% initialize S_int and S_ext
S_int = zeros( size( r ) );
S_ext = zeros( size( r ) );

% this loop evaluates Eq. (3.50) and (3.51)
for n = 0 : N - 1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N - 1 ) '.' ] );
   
    for m = -n : n
       
        % Eq. (2.38)
        S_breve = 4*pi .* 1i^( -abs( m ) ) .* sphharm( abs( m ), -m, phi_pw, theta_pw );
             
        % initialize G_breve_int and G_breve_ext
        G_breve_int_mm = 0;
        G_breve_int_nm = 0;
        G_breve_ext    = 0;
        
        % this loop evaluates Eq. (3.53) and (3.54)
        for n_prime = 0 : N_prime - 1
            for m_prime = -n_prime : n_prime
            
                G_breve_int_mm = G_breve_int_mm + G_breve_prime( 1 + n_prime^2 + n_prime + m_prime ) * (-1)^( abs( m ) + n_prime ) * ...
                        EI_sect( 1 + n_prime^2 + n_prime + m_prime, m + ( N - 1 ) + 1 );
                
                G_breve_int_nm = G_breve_int_nm + G_breve_prime( 1 + n_prime^2 + n_prime + m_prime ) * (-1)^( n + n_prime ) * ... 
                        EI( 1 + n^2 + n + m, 1 + n_prime^2 + n_prime + m_prime );
            
                G_breve_ext = G_breve_ext + G_breve_prime( 1 + n_prime^2 + n_prime + m_prime ) * (-1)^( n + n_prime ) * ...
                        EE( 1 + n^2 + n + m, 1 + n_prime^2 + n_prime + m_prime );
                    
            end
        end
       
        %%%%%%%%% TODO: skip rest of loop if G_breve_int_mm == 0 %%%%%%%%%%
        
        % Eq. (3.50) 
        S_int = S_int + S_breve .* G_breve_int_nm ./ G_breve_int_mm .* ...
                                sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );
   
        % Eq. (3.51)
        S_ext = S_ext + S_breve .* G_breve_ext ./ G_breve_int_mm .* ...
                             sphbesselh( n, 2, k.*r ) .* sphharm( n, m, beta, alpha );
                         
    end
    
end

toc;

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
