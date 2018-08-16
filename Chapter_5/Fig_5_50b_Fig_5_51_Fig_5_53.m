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

% Note that this script uses the myspectrogram function from 
% https://ccrma.stanford.edu/~jos/sasp/Matlab_listing_myspectrogram_m.html

clear;

figure_number = '5.50(b)'; % for Fig. 5.50(b)
%figure_number = '5.51(a)'; % for Fig. 5.51(a)
%figure_number = '5.51(b)'; % for Fig. 5.51(b)
%figure_number = '5.51(c)'; % for Fig. 5.51(c)
%figure_number = '5.51(d)'; % for Fig. 5.51(d)
%figure_number = '5.53(a)'; % for Fig. 5.53(a)
%figure_number = '5.53(b)'; % for Fig. 5.53(b)
%figure_number = '5.53(c)'; % for Fig. 5.53(c)
%figure_number = '5.53(d)'; % for Fig. 5.53(d)

if ( strcmp( figure_number, '5.53(a)' ) )
    f     = 1000;
    f_min = 700;
    f_max = 1300;
elseif ( strcmp( figure_number, '5.53(b)' ) )
    f     = 2000;
    f_min = 1400;
    f_max = 2600;
elseif ( strcmp( figure_number, '5.53(c)' ) )
    f     = 4000;
    f_min = 2800;
    f_max = 5200;
elseif ( strcmp( figure_number, '5.53(d)' ) )
    f     = 8000;   
    f_min = 5500;
    f_max = 10500; 
else
    f     = 500;
    f_min = 350;
    f_max = 650;
end

if ( ~strcmp( figure_number( 1 : 4 ), '5.51' ) )
    L = 50;
elseif ( strcmp( figure_number, '5.51(a)' ) )
    L = 10;
elseif ( strcmp( figure_number, '5.51(b)' ) )
    L = 20;    
elseif ( strcmp( figure_number, '5.51(c)' ) )
    L = 30;
elseif ( strcmp( figure_number, '5.51(d)' ) )
    L = 20;
end
 
if ( strcmp( figure_number, '5.51(d)' ) )
    tapering = true;
else 
    tapering = false;
end

v     = 40;
fs    = 44100;
c     = 343; 
omega = 2*pi*f;
M     = v / c;

dx_0  = 0.1;
d_ref = 1;

x = 0;
y = 4;
z = 0;

t = -1 : 1 / fs : 1;

% initialize s
s = zeros( size( t ) );

y_0 = 1;
z_0 = 0;

% loop over secondary sources (Eq. (5.69))
for x_0 = -L/2 : dx_0 : L/2

    disp( [ 'Processing secondary source at x_0 = ' num2str( x_0, '%0.1f' ) '.' ] ); 
    
    r_0 = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 + ( z - z_0 ).^2  );
    
    % Eq. (5.58) and (5.69)
    Delta = sqrt( ( x_0 - v * ( t - r_0/c ) ).^2  + ( y_0.^2 + z_0.^2 ) .* ( 1 - M^2 ) );

    % Eq. (5.57) and (5.69)
    tau = ( M .* ( x_0 - v * ( t - r_0/c ) ) + Delta ) / ( c * ( 1 - M^2 ) );

    % Eq. (5.66), (5.69), and (3.93)
    d = sqrt( 2*pi*d_ref / ( 1i * omega/c ) ) .* 2 * y_0 ./ Delta .* ( ( 1 - M^2 ) ./ Delta + 1/c * 1i*omega ) .* exp( 1i .* omega .* ( t - r_0/c - tau ) ) ./ Delta;

    if ( tapering )
        d = d * cos( x_0 / (L/2) * pi/2 )^2;
    end
    
    s = s + 1 / (4*pi) .* d ./ r_0;

end
    
% normalize
s = s ./ max( abs( s ) );

% align amplitude in the presence of artifacts
if ( strcmp( figure_number, '5.53(c)' ) )
    s = s .* 3;
elseif ( strcmp( figure_number, '5.53(d)' ) )
    s = s .* 4;
end

win_length = round( fs / 5 ); 
hop        = round( win_length / 32 ); 

spec      = myspectrogram( real( s ), win_length, fs, hamming( win_length ), win_length - hop, false );
size_spec = size( spec );

t_plot = linspace( t( 1 ), t( end ), size_spec( 2 ) );
f_plot = linspace(      0,       fs, size_spec( 1 ) );

figure;

imagesc( t_plot, f_plot, 20*log10( abs( spec ) / win_length ), [ -60 -20 ] );

turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
grid on;
axis square;

xlim(  [  -0.5   0.5 ] );
ylim(  [ f_min f_max ] );

xlabel( 't (s)' );
ylabel( 'f (Hz)' );

graph_defaults;
