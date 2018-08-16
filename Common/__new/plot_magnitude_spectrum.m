function [ ] = plot_magnitude_spectrum( spec, fs, new_figure )
%PLOT_MAGNITUDE_SPECTRUM plots lower half of magnitude spectrum on a 
% logarithmic scale
%   spec: complete spectrum
%   new_figure: true (default) or false

if ( nargin > 3 )
    new_figure = true;
end

spec = spec( 1 : end/2, : );

f = linspace( 5*eps, fs/2, size( spec, 1 ) );

if ( new_figure )
    figure;
else
    hold on;
end

semilogx( f, 20*log10( abs( spec ) ) );

if ( ~new_figure )
    hold off;
end

xlim( [ 30 fs ] );
grid on;

end

