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

subfigure = 'a'; % for Fig.4.14(a)
%subfigure = 'b'; % for Fig.4.14(b)

M = 80;
R = 1.5;

if ( subfigure == 'a' )
    r = R;
else
    r = R/2;
end

f    = linspace( 0, 3000, 200 ).';
k    = (2.*pi.*f)./343;
k(1) = k(2); % to avoid numerical instabilities

% initlialize G_ring_m
G_ring_m = zeros( length(f), 2*M+1 );

% this loop evaluates Eq. (2.34)
for m = -M : M
    
   disp( [ 'Calculating order ' num2str( m ) ' of ' num2str( M ) '.' ] );
    
    for n = abs(m) : M + 10 % go towards infinity

        G_ring_m( :, m+M+1 ) = G_ring_m( :, m+M+1 ) + ...
                    (-1i) .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, -m, pi/2, 0 ) .* ...
                            sphbesselj( n, k.*r ) .* sphharm( n, m, pi/2, 0 );
        
        % Note that setting m = 0 in the last sphharm function removes the
        % exponential function. This is easier than typing Eq. (2.34)
        % explicitely...
    end
    
end

figure;

if ( subfigure == 'a' )
    clims = [ -40 -20 ];
else
    clims = [ -80 -30 ];
end

imagesc( -M : M, f, 20*log10( abs( G_ring_m ) ), clims ); 

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
xlabel( 'm'      );
ylabel( 'f (Hz)' );
graph_defaults;
