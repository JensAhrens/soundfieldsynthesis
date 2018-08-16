function [ signal ] = convolve_signals( signal, varargin )
%CONVOLVE_SIGNALS Convolves SIGNAL with IMPULSE_RESPONSE
%   Uses fftfilt. Output SIGNAL has the same length like the input SIGNAL. 

for m = 1 : nargin-1
    
    impulse_response = varargin{ m };
    
    if ( size( impulse_response, 2 ) ~= 1 )
        error( 'Impulse response has to be single-channel.' );
    end
    
    % sometimes there are problems with short signals with fftfilt
    if ( size( signal, 2 ) < 2000 )
        for n = 1 : size( signal, 2 )
            signal( :, n ) = fftfilt( impulse_response, signal( :, n ) );
        end
    else
        for n = 1 : size( signal, 2 )
            signal( :, n ) = filter( impulse_response, signal( :, n ) );
        end
    end

end % for m = 1 : nargin-1

end

