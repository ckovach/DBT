

function [T,t,Err] = chopper(rg,evtt,fs) 

% [T,t,Err] = chopper(rg,evtt,fs) 
%
% Creates a matrix to segment a signal sampled at fs into windows ranging from rg(1) to
% rg(end) around event times specified in evtt.
%
% INPUT VARIABLES:
%
%  rg: 2 element array containing start and end times (e.g. -0.5 to 1.0) in seconds
%
%  evtt:  time stamps (in s)
%
%  fs:  sampling rate
%
%
% OUTPUT VARIABLE:
%
%  T: matrix of indices into a signal sampled at fs. Each column contains sample indices for 
%     the window surrounding the corresponding event time in evtt, as specified by rg. 
%     To segment the signal in the vector 'x' into a time-by-event matrix, use X = x(T), 
%     after correcting for any overlap with the ends of x. One way to implement such a 
%     correction is
%		T(T<1)=1; 
%		T(T>length(x)) = length(x);
%     which clamps values in X falling outside the time range of x to the first and last 
%     values of x, respectively.
%
%  t: vector, timestamps for the rows of T.
%  
%  Err: Difference between the sampled times and the exact window time.
% 

% C Kovach 2016

t = (rg(1):1/fs:rg(2))';

T = fs*( repmat(t,1,length(evtt)) + repmat(evtt(:)',length(t),1))+1;

if nargout > 2
   Err = (round(T)-T)./fs;  %%% Return the rounding error 
end

T = round(T);
