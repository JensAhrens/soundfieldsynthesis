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

% Note that the function Psi_extended_source.m used below uses the genHyper 
% function by Ben Barrowes. See file Common/genHyper_1.2/README for 
% credits and license.
%
% Note also that Fig. 5.23(a)-(c) are accidentally turned upside down in the 
% book. 
%
% Refer also to the comment on lines 42-44 below.

subfigure = 'a'; % for Fig. 5.23(a)
%subfigure = 'b'; % for Fig. 5.23(b)
%subfigure = 'c'; % for Fig. 5.23(c)

if ( subfigure == 'a' )
    eta       = 1;
    weight    = 1;
    alpha_rot = 0;
    
elseif ( subfigure == 'b' )
    eta       = 4;
    weight    = 1;
    alpha_rot = 0;
    
elseif ( subfigure == 'c' )  
    
    % Note that you can create your own complex source by choosing 
    % different eta, weight, and alpha_rot. They can be of arbitrary 
    % length, which has to equal for all three quantities though.
    
    eta       = [ 1 4 ];
    weight    = [ 1, .5 + 1.5*1i ]; 
    alpha_rot = [ 0 pi/8 ];
    
end    

N = 50;
a = 1;
f = 1000;
c = 343;
k = 2*pi*f/c;

% create spatial grid
X        = linspace( -3, 3, 300 );
Y        = linspace( -3, 3, 300 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% transfer to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x );
beta  = pi/2;

% initliaize S
S = zeros( size( x ) );

% this loop evaluates Eq. (5.27)
for n = 0 : N-1
    
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    hankel = sphbesselh( n, 2, k.*r ); 
    
    hankel_prime = 1 / ( 2*n + 1 ) * ( n * sphbesselh( n-1, 2, k*a ) - ( n + 1 ) * sphbesselh( n+1, 2, k*a ) ) * k;
    
    for m = -n : 2 : n
        
        % loop over all vibration modes
        for index = 1 : length( eta )
            
            % for n + abs(m) is odd
            if ( rem( n + abs( m ), 2 ) == 1 )
                V_breve = 0;

            else
                % Eq. (5.30)
                V_breve = (-1)^abs( m ) .* sqrt( ( 2*n + 1 ) / (4*pi) * ...
                        factorial( n - abs( m ) ) / factorial( n + abs( m ) ) ) * ...
                            Psi_extended_source( n, m ) .* Chi_extended_source( m, eta( index ) );            
            end

            % Eq. (5.27)
            S_breve = 1i .* V_breve ./ ( k .* hankel_prime );

            S = S + weight( index ) .* S_breve .* hankel .* sphharm( n, m, beta, alpha - alpha_rot( index ) );
            
        end % index
        
    end % m
    
end % n

S( r < a ) = 0;

% normalize
if ( subfigure == 'a' )
    norm_factor = 1000;
elseif ( subfigure == 'b' )
    norm_factor = 1000;
elseif ( subfigure == 'c' )
    norm_factor = 500;
end

S = S .* norm_factor;

figure;
imagesc( X, Y, real( S ), [ -1 1 ] );

% plot source
x_circ = linspace( -a, a, 100 );

hold on;
plot( x_circ,  sqrt( a.^2 - x_circ.^2 ), 'k', 'LineWidth', 2 )
plot( x_circ, -sqrt( a.^2 - x_circ.^2 ), 'k', 'LineWidth', 2 )
hold off;

turn_imagesc;
colormap gray;
axis square;
xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;
