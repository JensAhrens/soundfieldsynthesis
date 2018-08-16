function [ deviation, t ] = fit_line_to_ir_decay( t, decay, frame_size )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

deviation = zeros( size( decay ) );

%for n = 1 : frame_size : length( t ) - frame_size
for n = 1 : length( t ) - frame_size

    t_frame = t( n : n+frame_size-1 );
    decay_frame = decay( n : n+frame_size-1 );
    
    [ slope, offset ] = fit_line_to_data( t_frame, decay_frame );
    
    line = slope * t_frame + offset;
    
    %deviation( n : n+frame_size-1 ) = norm( decay_frame - line ) / sqrt( frame_size );
    deviation( n ) = norm( decay_frame - line ) / sqrt( frame_size );

end

