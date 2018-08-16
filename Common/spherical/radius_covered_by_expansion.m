function [ r_N ] = radius_covered_by_expansion( N, k )
% Calculates r_n according to Eq. (2.41)

r_N = N / k; 

if ( nargout == 0 )
    fprintf( 'r_N = %.3f\n', r_N );
end

end

