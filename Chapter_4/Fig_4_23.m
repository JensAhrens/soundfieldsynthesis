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

M_display = 80;
M         = 27;
R         = 1.5;

f    = linspace( 0, 3000, 200 ).';
k    = (2.*pi.*f)./343;
k(1) = k(2); % to avoid numerical instabilities

% initlialize G_ring_m
G_ring_m = zeros( length(f), 2*M_display+1 );

% this loop evaluates a part of Eq. (2.34)
for m = -M_display : M_display
    
   disp( [ 'Calculating order ' num2str( m ) ' of ' num2str( M_display ) '.' ] );
    
    for n = abs(m) : M_display + 10 % go towards infinity

        G_ring_m( :, m+M_display+1 ) = G_ring_m( :, m+M_display+1 ) + ...
                    (-1i) .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, -m, pi/2, 0 ) .* ...
                            sphbesselj( n, k.*R ) .* sphharm( n, m, pi/2, 0 );
        
        % Note that setting m = 0 in the last sphharm function removes the
        % exponential function. This is easier than typing Eq. (2.34)
        % explicitely...
    end
    
end

% Fig. 4.23(a)
figure;

imagesc( -M_display : M_display, f, 20*log10( abs( G_ring_m ) ), [ -40 -20 ] ); 

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
xlabel( 'm'      );
ylabel( 'f (Hz)' );
graph_defaults;

G_ring_m( :, 1 : -M + M_display )      = 0;
G_ring_m( :, M + M_display + 2 : end ) = 0;

% Fig. 4.23(b)
figure;

imagesc( -M_display : M_display, f, 20*log10( abs( G_ring_m ) ), [ -40 -20 ] ); 

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
xlabel( 'm'      );
ylabel( 'f (Hz)' );
graph_defaults;
