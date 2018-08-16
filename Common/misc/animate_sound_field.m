function [ ] = animate_sound_field( X, Y, S, clims, T, t_wait )
%ANIMATE_SOUND_FIELD Summary of this function goes here
%   Detailed explanation goes here

t = 1;

figure;
while ( 1 )

    imagesc( X, Y, real( S .* exp( 1i .* t/T .* 2*pi ) ) );

    colormap jet;
    turn_imagesc;
    axis square;
    axis equal;
    axis tight;
    caxis( clims );
    drawnow;
    
    pause( t_wait );
    
    t = t + 1;
end

