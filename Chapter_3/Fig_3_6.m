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

% both Fig. 3.6(a) and (b) will be created 

N = 3;

alpha_s = 0;
beta_s  = pi/2;

% create spatial grid
Alpha = linspace( -pi, pi, 201 );
Beta  = linspace(   0, pi, 101 );

[ alpha, beta ] = meshgrid( Alpha, Beta );

D = zeros( size( alpha ) );

% this loop evaluates Eq. (3.35)
for n = 0 : N-1
    for m = -n : n

        D = D + sqrt( 2*n+1 ) .* sphharm( n, -m, beta_s, alpha_s ) ./ ...
                                 sphharm( n, 0, 0, 0 ) .* ...
                                 sphharm( n, m, beta, alpha );
        
    end
end

% suppress imaginary component that arise due to limited numeric precision
D = real( D );

% normalize
D = D ./ max( abs( D( : ) ) );

%figure;
imagesc( Alpha / pi, Beta / pi, D, [ -.5 1 ] );
colormap gray;
revert_colormap;
colorbar;
xlabel('$\alpha / \pi$', 'Interpreter', 'latex' );
ylabel('$\beta  / \pi$', 'Interpreter', 'latex' );
graph_defaults;

figure;
plot( Alpha / pi, real( D( (end+1)/2, : ) ), 'Color', [ .5 .5 .5 ], 'Linewidth', 2 );
ylim( [ -.5 1 ] );
grid on;
xlabel('$\alpha / \pi$', 'Interpreter', 'latex' );
graph_defaults;
