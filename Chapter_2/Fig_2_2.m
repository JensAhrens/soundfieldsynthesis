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

subfigure = 'a'; % for Fig. 2.2(a)
%subfigure = 'b'; % for Fig. 2.2(b)
%subfigure = 'c'; % for Fig. 2.2(c)

N = 5; % order
x = linspace( 0, 10, 1000 ); % argument of functions

% intialize output variable
out = zeros( N+1, length( x ) );


for n = 0 : N
    if ( subfigure == 'a' )
        out( n+1, : ) = sphbesselj( n, x );
    elseif ( subfigure == 'b' )
        out( n+1, : ) = sphbessely( n, x );
    elseif ( subfigure == 'c' )
        out( n+1, : ) = 20*log10( abs( sphbesselh( n, 2, x ) ) );
    end
end

% create different color shadings for plots
color_order = repmat( linspace( 0, 1, N+2 ).', [ 1 3 ]);
    
figure;
set( gca, 'ColorOrder', color_order      );
set( gca, 'NextPlot'  , 'replacechildren');

plot( x, out, 'LineWidth', 2 );

if ( subfigure == 'a' )
    axis( [ x( 1 ) x( end ) -0.5 1.5 ] );
elseif ( subfigure == 'b' )
    axis( [ x( 1 ) x( end ) -100 0 ] );
elseif ( subfigure == 'c' )
    axis( [ x( 1 ) x( end ) -50 100 ] );
end

grid on;
xlabel( 'kr' )
graph_defaults;
