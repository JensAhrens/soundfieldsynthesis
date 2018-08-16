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

subfigure = 'a'; % for Fig. 4.6(a)
%subfigure = 'b'; % for Fig. 4.6(b)

N    = 81;
f    = linspace( 0, 3000, 200 ).';
f(1) = f(2); % to avoid numerical instabilities
k    = (2.*pi.*f)./343;
R    = 1.5;

if ( subfigure == 'a' )
    r = 3/4*R;
else
    r = R/2;
end

% initialize G_ring_n0
G_ring_n0 = zeros( length(k), N );

for n = 0 : N - 1
    % from Eq. (2.37a)
    G_ring_n0( :, n+1 ) = (-1i) .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, 0, 0, 0 ) .* sphbesselj( n, k.*r );         

end


figure;
imagesc( 0 : N - 1, f, 20*log10( abs( G_ring_n0 ) ), [ -70 -10 ] );
turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
xlabel( 'n'      );
ylabel( 'f (Hz)' );
graph_defaults;
