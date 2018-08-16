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

% both Fig. 5.7(a) and (b) will be created 

r_s     = 3;
alpha_s = -pi/2;

N       = 40;
R       = 1.5;  % m
f       = 1000; % Hz
c       = 343;  % m/s

k = (2*pi*f)/c;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x );
beta  = pi/2;

% initialize S_int and S_ext
S_int = zeros( size( r ) );
S_ext = zeros( size( r ) );

% this loop evaluates Eq. (3.50) and (3.51)
for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
       
    for m = -n : n
        
        % Eq. (2.37a)
        S_breve = -1i .* k .* sphbesselh( abs(m), 2, k.*r_s ) .* sphharm( abs(m), -m, pi/2, alpha_s );   
        
        % Eq. (2.37a)
        G_breve_int_nm = -1i .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, -m, pi/2, 0 );

        % Eq. (2.37a)
        G_breve_int_mm = -1i .* k .* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 );

        % Eq. (2.37b)
        G_breve_ext    = -1i .* k .* sphbesselj( n, k.*R ) .* sphharm( n, -m, pi/2, 0 );
    
        % Eq. (3.50)
        S_int = S_int + S_breve .* G_breve_int_nm ./ G_breve_int_mm .* ...
                                sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );
        
        % Eq. (3.51)
        S_ext = S_ext + S_breve .* G_breve_ext ./ G_breve_int_mm .* ...
                                sphbesselh( n, 2, k.*r ) .* sphharm( n, m, beta, alpha );      
        
     end
end

% assemble exterior and interior field
S          = zeros( size( x ) ); 
S( r < R ) = S_int( r < R );
S( r > R ) = S_ext( r > R );

% normalize
S = S ./ abs( S( end/2, end/2 ) );

% Fig. 5.7(a)
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

% Fig. 5.7(b)
figure;
imagesc( X, Y, 20*log10 ( abs( S ) ), [ -10 10 ] );
turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;

hold on;

% plot secondary source distribution
plot( x_circ,  sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )
plot( x_circ, -sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;
