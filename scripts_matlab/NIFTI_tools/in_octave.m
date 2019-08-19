%
% Simple function to check if we are are currently executing in Octave or Matlab   
%
% USAGES: inOctave = in_octave();
%         [inOctave,version_str] = in_octave();
%         
%  
% Mark Prati (LMP)- mprati@research.baycrest.org 
%
%
function [inOctave,version_str] = in_octave()

try
    ver_num = OCTAVE_VERSION;
    inOctave = 1;
    version_str = ['OCTAVE ' ver_num];
catch
    inOctave = 0;
    version_str  = ['MATLAB ' version];
end
