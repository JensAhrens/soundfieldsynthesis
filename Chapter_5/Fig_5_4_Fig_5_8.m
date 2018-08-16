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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The author of this script is Sascha Spors.

clear;

subfigure = '5.4'; % for Fig. 5.4
%subfigure = '5.8'; % for Fig. 5.8

ns   = [ 1 5 10 20 28 ]; % orders to plot
r_s  = 3; % for Fig. 5.8
R    = 1.5;
c    = 343;
fs   = 44100;
f    = linspace( 0, fs/2 - 1/fs, 1024 );
f(1) = f(2); % avoid numeric instabilities
k    = 2*pi*f/c;


for n = 1 : length( ns )

    % analytic reference
    if ( strcmp( subfigure, '5.4' ) )
        
        h1( n, : ) = 2 * 1i * 1i^( -ns(n) ) ./ ( k .* R .* sphbesselh( ns(n), 2, k*R ) );
        src_type   = 'pw';
    
    elseif ( strcmp( subfigure, '5.8' ) )
        
        h1( n, : ) = sphbesselh( ns(n), 2, k*r_s ) ./ sphbesselh( ns(n), 2, k*R ); 
        src_type   = 'sw';
        
    end
        
    %  modal filter design
    hd = HOA25D_modal_filter( src_type, R, ns(n), fs, c, r_s );
    
    % compute frequency response of resulting filter
    [ h2( n, : ), f2 ] = freqz( hd, length( f ), fs );
    h2( n, : ) = h2( n, : ) .* exp( 1i * 2 * pi * f * ( r_s - R ) / c );
    
end
   
if ( strcmp( subfigure, '5.4' ) )
    % avoid log of 0
    h2 = h2 + 5*eps;
elseif ( strcmp( subfigure, '5.8' ) )
    % manual normalization
    h2 = h2 .* 4.7;
end

figure;

ph1 = semilogx( f, 20*log10( abs ( h1( 1, : ).' ) ), 'Color', [ 0.1667    0.1667    0.1667 ], 'LineWidth', 2 );
hold on;
      semilogx( f, 20*log10( abs ( h2( 1, : ).' ) ), 'k:', 'LineWidth', 2 );
ph2 = semilogx( f, 20*log10( abs ( h1( 2, : ).' ) ), 'Color', [ 0.3333    0.3333    0.3333 ], 'LineWidth', 2 );
      semilogx( f, 20*log10( abs ( h2( 2, : ).' ) ), 'k:', 'LineWidth', 2 );
ph3 = semilogx( f, 20*log10( abs ( h1( 3, : ).' ) ), 'Color', [ 0.5000    0.5000    0.5000 ], 'LineWidth', 2 );
      semilogx( f, 20*log10( abs ( h2( 3, : ).' ) ), 'k:', 'LineWidth', 2 );
ph4 = semilogx( f, 20*log10( abs ( h1( 4, : ).' ) ), 'Color', [ 0.6667    0.6667    0.6667 ], 'LineWidth', 2 );
      semilogx( f, 20*log10( abs ( h2( 4, : ).' ) ), 'k:', 'LineWidth', 2 );
ph5 = semilogx( f, 20*log10( abs ( h1( 5, : ).' ) ), 'Color', [ 0.8333    0.8333    0.8333 ], 'LineWidth', 2 );
      semilogx( f, 20*log10( abs ( h2( 5, : ).' ) ), 'k:', 'LineWidth', 2 );

hold off;

axis( [ 20 20000 -150 50 ] );

xlabel( 'f (Hz)' );
ylabel( 'mag'    );
grid on;

phl = legend( [ ph1(1) ph2(1) ph3(1) ph4(1) ph5(1) ], '$m$ = 1', '$m$ = 5', '$m$ = 10', '$m$ = 20', '$m$ = 28', 'Location', 'SouthEast' );
set( phl, 'Interpreter', 'latex' )
graph_defaults;
