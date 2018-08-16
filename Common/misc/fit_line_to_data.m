function [ slope, offset ] = fit_line_to_data( x, y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = x(:);
x = [ x, ones( size( x ) ) ];

y = y(:);

% Ax = b
param = x \ y;

slope  = param( 1 );
offset = param( 2 );

end

