function loss = ideloss(ry,t,xn) 
% input
%   ry: r and eta_y
%   t:  time instants
%   xn: observations at time instants of t

dt = diff(t);
n = length(xn);

wn = cumsum(dt.*(xn(1:end-1)+xn(2:end)))/2;
Theta = [ry(2) + wn, (ry(2) + wn).^ry(1)];

loss  = (eye(n-1) - Theta*pinv(Theta) )*xn(2:n);

assignin('base','a12', Theta\xn(2:n) );

end
