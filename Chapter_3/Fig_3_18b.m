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

N        = 3; % order of directivity; the books states N = 13 but the simulation actually shows N = 3 (sorry...)
y_ref    = 1;
theta_pw = pi/4;
phi_pw   = pi/2;
f        = 1000; % Hz
c        = 343;  % m/s
k        = (2*pi*f)/c;

k_pw_x = k * cos( theta_pw ) * sin( phi_pw );
k_pw_y = k * sin( theta_pw ) * sin( phi_pw );
k_pw_z = k * cos( phi_pw );

%%%%%%%%%%%%%%%%%%%%%%%%%% prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_interval = [ -50 50 ];
delta_x          = .01; % sampling interval for spatial fft in meters

X = spatial_interval( 1 ) : delta_x : spatial_interval( 2 );
Y = linspace( -1, 3, 201 );

k_x_s = (2*pi) / delta_x; % spatial sampling frequency

% create k_x
k_x    = linspace( 0, k_x_s/2, ( length( X ) + 1 ) / 2 ); % positive frequencies
k_x(1) = k_x(2); % to avoid numerical instabilities
k_x    = [ -fliplr( k_x( 2 : end ) ), k_x ]; % adds negative frequencies

% create 2D grid
[ k_x_m, y_m ]    = meshgrid( k_x, Y );
y_m               = abs( y_m );
y_m( y_m < 0.01 ) = 0.01; % to avoid numerical instablilities
%%%%%%%%%%%%%%%%%%%%%%% end prepare spatial fft %%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% calculate G_tilde %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = atan2( y_ref, X );
r     = sqrt( X.^2 + y_ref.^2 );
G     = zeros( size( r ) );

% this loop evaluates Eq. (2.32b) with Eq. (2.44) along the reference line
for n = 0: N-1
    spherical_hankel = sphbesselh( n, 2, k.*r );
    
    for m = -n : n
        G = G + i^(-n) .* ( factorial( N-1 ) * factorial( N ) ) / ( factorial( N+n ) * factorial( N-n-1 ) ) .* ...
                                 sphharm( n, -m, pi/2, pi/2 ) .* spherical_hankel .* sphharm( n, m, pi/2, alpha );
    end
    
end

G_tilde = fftx( G, [], 2 );

% determine index of k_pw_x in k_x 
[ tmp ind ] = find( k_x < k_pw_x );
index = ind( end );

% this mimicks 2pi dirac( k_x - k_pw_x ) in Eq. (C.5)
G_tilde_yref = 2*pi .* G_tilde( index );

% free some memory
clear alpha r G G_tilde;
%%%%%%%%%%%%%%%%%%%% end calculate G_tilde %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% calculate synthesized sound field %%%%%%%%%%%%%%%%%%%%%%
% spatial grid
[ x, y ] = meshgrid( X, Y );
alpha    = atan2( y, x );
r        = sqrt( x.^2 + y.^2 );

G = zeros( size( r ) );

% calculate G; this loop evaluates Eq. (2.32b) with Eq. (2.44)
for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    spherical_hankel = sphbesselh( n, 2, k.*r );
    
    for m = -n : n
        G = G + 1i^(-n) * ( factorial( N-1 ) * factorial( N ) ) / ( factorial( N+n ) * factorial( N-n-1 ) ) .* ...
                 sphharm( n, -m, pi/2, pi/2 ) .* spherical_hankel .* sphharm( n, m, pi/2, alpha );
    end
    
end

% spatial Fourier transform
G_tilde = fftx( G, [], 2 );

% this is actually an inverse (time) Fourier transform of Eq. (3.71) with
% (3.77)
S = exp( -1i .* k_pw_y .* y_ref ) .* exp( -1i .* k_pw_x .* x ) .* ...
           repmat( G_tilde( :, index ), [ 1 size(G_tilde, 2) ] ) ./ G_tilde_yref;

% normalization
S = S ./ abs( S( (end+1)/2, (end+1)/2 ) ); 

figure;
imagesc( X, Y, real( S ), [ -2 2 ] )
colormap gray;
xlim( [ -2 2 ] );
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

hold on;
% plot secondary source distribution
plot( [ -2 2 ], [ 0 0 ], 'k', 'LineWidth', 2 );
hold off;

graph_defaults;
