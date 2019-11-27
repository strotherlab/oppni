function [errmsg, msgobj] = freqchk_octave(wp,ws,opt)
% Modified version of freqcheck from the MATLAB stats toolbox, so that it
% works with the default settings in octave. The defaults must be used for
% now.
% Edits By Andrew Lofts July 4 2018
%
%
%FREQCHK Parameter checking for BUTTORD, CHEB1ORD, CHEB2ORD, ELLIPORD.
%
%   MSG=FREQCHK(WP,WS,OPT) checks for the correct passband (WP) and 
%   stopband (WS) frequency specifications and returns a diagnostic 
%   message (MSG). If the parameters are specified correctly MSG is
%   empty. OPT is a string indicating digital or analog filter design.
%
%   The frequency parameters are checked to be of the same length.
%   For bandpass and bandstop filters, the frequency bands are checked
%   for increasing values and no overlaps.
% 
%   For digital filters, WP and WS are checked to be in (0,1). 
%   For analog filters, WP and WS are checked to be non-negative values.
   
%   Copyright 1988-2002 The MathWorks, Inc.

errmsg='';
msgobj=[];

% Check for correct lengths
if length(wp)~=length(ws),
   msgobj = error('Missmatched Frequency Lengths, see MATLAB version');
   errmsg = getString(msgobj);
   return
end

% Check for allowed interval values 
if strcmp(opt,'z'),
   if any(wp<=0) || any(wp>=1) || any(ws<=0) || any(ws>=1),
      msgobj = error('Cutoff outside of the frequency see MATLAB version');
      errmsg = getString(msgobj);
      return    
   end
else % Analog filter design
   if any(wp<=0) || any(ws<=0),
      msgobj = error('Analog error see MATLAB version');
      errmsg = getString(msgobj);
      return    
   end
end

% For Band specifications
if length(wp)==2,
   % Check for frequencies to be in increasing order
   if (wp(1)>=wp(2)) || (ws(1)>=ws(2)),
      msgobj = error('Cutoff issue see MATLAB version');
      errmsg = getString(msgobj);
      return
   end    
   % Check for passband and stopband frequency overlaps  
   if ~( ((wp(1)<ws(1)) && (wp(2)>ws(2))) || ((wp(1)>ws(1)) && (wp(2)<ws(2))) ),
      msgobj = error('No Overlaps see MATLAB version');
      errmsg = getString(msgobj); 
      return
   end
end

% [EOF] freqchk.m

