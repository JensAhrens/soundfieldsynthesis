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

subfigure = 'a'; % for Fig. 2.6(a)
%subfigure = 'b'; % for Fig. 2.6(b)

N     = 81; % summation limit 
f_max = 3000;

f     = linspace( 0, f_max, 200 ).';
f(1)  = f(2); % to avoid numerical problems
c     = 343;
k     = ( 2.*pi.*f ) ./ c;
R     = 1.5;
r     = R/2;

S_ring_n = zeros( length( f ), N );

for n = 0 : N-1
    
    if ( subfigure == 'a' )
        % plane wave
        S_ring_n( :, n+1 ) = 4*pi*( 1i )^(-n) .* ...
                    sphharm( n, 0, pi/2, 0 ) .* sphbesselj( n, k.*r );
    elseif ( subfigure == 'b' )
        % point source
        S_ring_n( :, n+1 ) = ( -1i ) .* k .* sphbesselh( n, 2, k.*R ) .* ...
                    sphharm( n, 0, pi/2, 0 ) .* sphbesselj( n, k.*r );         
    end
        
end % n


figure;
imagesc( ( 0 : N-1 ), f, 20*log10( abs( S_ring_n ) ) ); 

if ( subfigure == 'a' )
    caxis( [ -40 0 ] );
elseif ( subfigure == 'b' )
    caxis( [ -60 -20 ] );
end

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;

hold on;
% plot black line along f = ( n*c ) / ( 2*pi*r )
plot( ( 0 : N-1 ), ( 0 : N-1 ) * c / ( 2*pi*r ), 'k', 'Linewidth', 2, 'Color', [ 0 0 0 ] );
hold off;

axis square;
xlabel( 'n' );
ylabel( 'f (Hz)' );
graph_defaults;
