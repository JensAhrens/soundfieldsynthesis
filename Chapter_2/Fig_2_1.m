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

subfigure = 'a'; % for Fig. 2.1 (a)
%subfigure = 'b'; % for Fig. 2.1 (b)

% frequency of plane wave
f = 1000;
c = 343;
k = 2*pi*f/c;

if ( subfigure == 'a' )
    k_pw = [ k 0 0 ];                        
    
elseif ( subfigure == 'b' )
    k_pw = [ k .* sqrt( 1.01 ) -1i.*k .* sqrt( .01 ) 0 ]; 
    
end

% create a spatial grid
resolution = 300;
X          = linspace( -2, 2, resolution );
Y          = linspace(  0, 4, resolution );
[ x, y ]   = meshgrid( X, Y );
z          = 0;

% Eq. (2.7)
S = exp( -1i .* k_pw( 1 ) .* x ) .* exp( -1i .* k_pw( 2 ) .* y ) .* exp( -1i .* k_pw( 3 ) .* z );

% normalize
S = S ./ abs( S( end/2, end/2 ) );

figure;
imagesc( X, Y, real( S ) ); 

if ( subfigure == 'a' )
    caxis( [ -2  2 ] );
elseif ( subfigure == 'b' )
    caxis( [ -20 20 ] );
end

turn_imagesc;
colormap gray;
axis square;
xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;


