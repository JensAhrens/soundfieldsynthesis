function [ imp_resp, spec ] = ls_system_id( reference_signal, measured_signal, imp_resp_length )
% LS_SYSTEM_ID performs linear system ID from reference signal to measured
%              signal using virtual overlap-save convolution

if ( nargin > 3 )
    imp_resp_length = 1024;
end

if ( size( reference_signal, 2 ) ~= 1 || size( measured_signal, 2 ) ~= 1 )
    error( 'Only single-channel signals are possible.' );
end

% assure proper length of signals
if ( size( reference_signal, 1 ) > size( measured_signal, 1 ) )
	measured_signal  = [ measured_signal; zeros( size( reference_signal, 1 ) - size( measured_signal, 1 ), 1 ) ]; 
else
    reference_signal = [ reference_signal; zeros( size( measured_signal, 1 ) - size( reference_signal, 1 ), 1 ) ]; 
end

reference_signal = [ zeros( imp_resp_length, 1 ); reference_signal ];
measured_signal  = [ measured_signal; zeros( imp_resp_length, 1 ) ];

%%%%%%%%%%%%%% divide signal into frames and do FFT %%%%%%%%%%%%%%%%%%%%%%%

% allocate memory
reference_signal_spectrum = zeros( 2*imp_resp_length, 0 );
measured_signal_spectrum  = zeros( 2*imp_resp_length, 0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1 : imp_resp_length : size( reference_signal, 1 ) - 2*imp_resp_length + 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    reference_signal_spectrum = [ reference_signal_spectrum, fft( reference_signal( n : n + 2*imp_resp_length - 1, 1 ), [], 1 ) ];
    measured_signal_spectrum  = [ measured_signal_spectrum , fft( measured_signal(  n : n + 2*imp_resp_length - 1, 1 ), [], 1 ) ];
    
end

%%%%%%%%%%%%%%%%%%%%%% determine LS solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning( 'Not implemented correctly.' );
spec = reference_signal_spectrum.' \ measured_signal_spectrum;

imp_resp = ifft( spec, [], 1 );

%%%%%%%%%%%%%%%%%%%% remove virtual zero-padding %%%%%%%%%%%%%%%%%%%%%%%%%%
tail     = imp_resp( end/2+1 : end, 1 );
imp_resp = imp_resp( 1 : end/2, 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%% sanity check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMS_imp  = norm( imp_resp ) / sqrt( length( imp_resp ) );
RMS_tail = norm( tail     ) / sqrt( length( tail     ) );

fprintf( '\nLS system ID:\n' );
fprintf( 'RMS imp resp: %0.2f\n', RMS_imp );
fprintf( 'RMS tail    : %0.2f (should tend to 0)\n\n', RMS_tail );

spec = fft( imp_resp, [], 1 );

end

