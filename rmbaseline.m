
function out = rmbaseline(bx,usepts)


% Fits a polynomial to the spectrum and normalizes by the fitted value.

% C Kovach 2013
% 
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 2 || isempty(usepts)
    usepts = true;
end
polyord = 8;

mn = mean(abs(bx.blrep));
kp = find(mn~=0 & usepts);
p = polyfit(kp./length(mn),log(mn(kp)),polyord);

bl = polyval(p,(1:length(mn))/length(mn));



out = bx.blrep./repmat(exp(bl),size(bx.blrep,1),1);

