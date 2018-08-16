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

subfigure = 'a'; % for Fig. 4.20(a)
%subfigure = 'c'; % for Fig. 4.20(c)

theta_pw = pi/2;
phi_pw   = pi/2;
L        = 56;
M        = 27;
R        = 1.5;
x        = 1;
y        = 0;

no_of_frequency_bins = 1025;

fs   = 44100;
c    = 343;
f    = linspace( 0, fs, no_of_frequency_bins + 1 );
f    = f( 1 : end/2 ).';
k    = 2.*pi.*f./c;
k(1) = k(2); % to avoid numerical instabilities

normalization_factor = 0.7; % amplitude normalization
cutoff               = 2200; % for Fig. 4.20(c)

% time instances to plot
t_to_show    = [ -0.004 0.01 ]; % in sec
taps_to_show = round( ( t_to_show(1) * fs : t_to_show(2) * fs ) + 1.5 * no_of_frequency_bins );

d_alpha_0 = 2*pi/L;

% initialize S_breve and G_breve
S_breve = zeros( length(f), 2*M+1 );
G_breve = zeros( length(f), 2*M+1 );

for m = -M : M

    % Eq. (2.38)
    S_breve( :, m+M+1 ) = 1i.^( -abs(m) ) .* sphharm( abs(m), -m, phi_pw, theta_pw );
    
    % Eq. (2.37a)
    G_breve( :, m+M+1 ) = -1i .* k .* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 );
    
end

% initialize D
D = zeros( length(f), L );

% loop over all secondary sources
for l = 0 : L-1
    
    % this loop evaluates Eq. (3.49)
    for m = -M : M    
        alpha_0    = l * d_alpha_0;
        D( :,l+1 ) = D( :, l+1 ) + 1 / ( 2*pi*R ) .* ...
                S_breve( :, m+M+1 ) ./ G_breve( :, m+M+1 ) .* exp( i .* m .* alpha_0 );    
        
    end
end

d = real( ifft( [ D; conj( flipud( D( 2 : end-1, : ) )  ) ] ) );
% move t = 0 to the center of the buffer
d = circshift( d, [ length(f)-1 0 ] );

% put zeros around to have some headroom
d = [ zeros( no_of_frequency_bins, L ); d; zeros( no_of_frequency_bins, L ) ];

% normalize
d = d ./ max(abs( d( : ) ) );

% initialize s
s = zeros( length( taps_to_show ), 1 );

% this loop evaluates Eq. (5.69)
for l = 1 : L
    
    x_0 = R * cos( (l-1) * d_alpha_0 );
    y_0 = R * sin( (l-1) * d_alpha_0 );
    
    r = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 );
    
    % from Eq. (5.69)
    t = ( taps_to_show./fs - r./c ) .* fs + 1; % in samples

    % Interpolate the impulse responses to find the values at instances
    % t that we are interested in. 
    d_interp = interp1( ( 1 : size( d, 1 ) ), d( :, l ), t, 'linear').';
    
    % from Eq. (5.69)
    s = s + d_interp ./ r;
    
end

% normalize
s = s ./ normalization_factor;

figure;

if ( subfigure == 'a' )
    % avoid log of 0
    s = 20*log10( abs( s ) + 5*eps );
    
    plot( ( taps_to_show - 1.5 * no_of_frequency_bins ) ./ fs .* 1000, s, 'k' )
    
else
    % apply filtering in FFT domain (ugly but serves its purpose...)
    S   = fft( s );
    S_1 = S;
    S_2 = S;
    
    % cutoff frequency in samples
    cutoff = round( cutoff./fs .* size( s, 1 ) ) + 1; % samples

    % set all portions that we want to suppress to zero
    S_1( 1 : cutoff )           = 0; 
    S_1( end-cutoff+1 : end )   = 0;
    S_2( cutoff+1 : end-cutoff) = 0;

    s_1 = real( ifft( S_1 ) );
    s_2 = real( ifft( S_2 ) );
    
    % avoid log of 0
    s_1 = 20*log10( abs( s_1 ) + 5*eps );
    s_2 = 20*log10( abs( s_2 ) + 5*eps );

    plot( ( taps_to_show - 1.5 * no_of_frequency_bins ) ./ fs  .* 1000, s_1, 'k', 'Color', [ .7 .7 .7 ] )
    hold on;
    plot( ( taps_to_show - 1.5 * no_of_frequency_bins ) ./ fs  .* 1000, s_2, 'k', 'LineWidth', 2 )
    hold off;

    legend( 'hp', 'lp', 'Location', 'NorthEast' );
    
end
    
ylim( [ -30 0 ] );
xlim( t_to_show .* 1000 );
axis square;
grid on;
xlabel( 't (ms)' );
graph_defaults;
