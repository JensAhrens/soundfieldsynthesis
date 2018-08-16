function [ signal, offset ] = apply_fractional_delay( signal, fractional_delay, method, N )
% Applies the fractional delay fractional_delay ( < 1 ) to single-channel 
% signal; output has same length like input, may be mutlichannel
% fractional_delay: may be in interval [0, 1]
% method: 'lagrange' (default)
% N     : interpolation order, default: 5
%
% offset: inherent integer delay applied by fractional procedure
%
% Requires FDTOOLS package by Laakso, Välimäki, Karjalainen, and Laine

if ( nargin < 3 )
    method = 'lagrange';
end
if ( nargin < 4 )
    N = 5;
end

% include suitable integer delay
% N odd
if ( rem( N, 2 ) )
    offset = ( N - 1 ) / 2;
% N even
else
    offset = N / 2;
end

if ( fractional_delay > N )
    error( 'Fractional delay has to be smaller than N' );
end

switch method
    case 'lagrange'
        h = lagrange( N, fractional_delay );
        
        % loop over channels
        for n = 1 : size( signal, 2 )
            signal( :, n ) = filter( h, 1, signal( :, n ) );
        end
        
    otherwise
        error( 'Unknown method.' )
        
end


end

function h = lagrange( N, delay )
%LAGRANGE  h=lagrange(N,delay) returns order N FIR 
%          filter h which implements given delay 
%          (in samples).  For best results, 
%          delay should be near N/2 +/- 1.
% 
% Downloaded from https://ccrma.stanford.edu/~jos/Interpolation/
% Matlab_Code_Lagrange_Fractional.html

n = 0:N;
h = ones(1,N+1);
for k = 0:N
    index = find(n ~= k);
    h(index) = h(index) *  (delay-k)./ (n(index)-k);
end

end

