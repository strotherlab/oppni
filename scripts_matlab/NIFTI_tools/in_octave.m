%
% Simple function to check if we are are currently executing in Octave or Matlab   
%
% Mark Prati - mprati@research.baycrest.org 
%
function inOctave = in_octave()

try
    ver_num = OCTAVE_VERSION;
    inOctave = 1;
    version_str = ['OCTAVE ' ver_num];
catch
    inOctave = 0;
    version_str  = ['MATLAB ' version];
end
