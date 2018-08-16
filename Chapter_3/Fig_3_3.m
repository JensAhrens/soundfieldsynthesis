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

% both Fig. 3.3(a) and (b) will be created 

theta_pw = pi/2;
phi_pw   = pi/2;
N        = 40;
R        = 1.5;  % m
f        = 1000; % Hz
c        = 343;  % m/s

k = (2*pi*f)/c;

% create spatial grid
Y        = linspace( -2, 2, 200 );
Z        = linspace( -2, 2, 200 );
[ y, z ] = meshgrid( Y, Z );
x        = 0;

% convert to spherical coordinates
r     = sqrt( y.^2 + z.^2 );
alpha = atan2( y, x );
beta  = acos( z ./ r );

% initialize S_int and S_ext
S_int = zeros( size( r ) );
S_ext = zeros( size( r ) );

% this loop evaluates Eq. (3.22) and (3.23)
for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
   
    for m = -n : n
       
        % Eq. (2.38)
        S_breve = 4*pi .* 1i^( -n ) .* sphharm( n, -m, phi_pw, theta_pw );
       
        S_int = S_int + S_breve .* sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );

        % from Eq. (2.37a)
        G_breve_int = -1i .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, 0, 0, 0 );
       
        % from Eq. (2.37b)
        G_breve_ext = -1i .* k .* sphbesselj( n, k.*R )    .* sphharm( n, 0, 0, 0 );
        
        S_ext = S_ext + S_breve .* G_breve_ext ./ G_breve_int .* ...
                             sphbesselh( n, 2, k.*r ) .* sphharm( n, m, beta, alpha );
        
    end
    
end

% assemble exterior and interior field
S          = zeros( size( r ) ); 
S( r < R ) = S_int( r < R );
S( r > R ) = S_ext( r > R );

figure;
imagesc( Y, Z, real( S ), [ -2 2 ] );
turn_imagesc;
colormap gray;
axis square;

hold on;

% plot secondary source distribution
y_circ = linspace( -R, R, 200 );
plot( y_circ,  sqrt( R.^2 - y_circ.^2 ), 'k', 'Linewidth', 2 )
plot( y_circ, -sqrt( R.^2 - y_circ.^2 ), 'k', 'Linewidth', 2 )

% plot horizontal plane
plot( [ -2 2 ], [ 0 0 ], 'k--' );

hold off;

xlabel( 'y (m)' );
ylabel( 'z (m)' );
graph_defaults;

% this avoids wigglyness due to numercial artifacts 
%S( r < R ) = 1;

figure;
imagesc( Y, Z, 20*log10 ( abs( S ) ), [ -10 10 ] );
turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;

hold on;

% plot secondary source distribution
plot( y_circ,  sqrt( R.^2 - y_circ.^2 ), 'k', 'Linewidth', 2 )
plot( y_circ, -sqrt( R.^2 - y_circ.^2 ), 'k', 'Linewidth', 2 )

% plot horizontal plane
plot( [ -2 2 ], [ 0 0 ], 'k--' );

hold off;

xlabel( 'y (m)' );
ylabel( 'z (m)' );
graph_defaults;
