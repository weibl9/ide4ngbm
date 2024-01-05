function loss = odeloss(ry,t,yn)
% input
%   ry: r and eta_y
%   t:  time instants
%   yn: cusum observations [eta_y; eta_y+cusum]

n = length(yn);
zn = diff(yn)./diff(t);

Theta = [yn(1:end-1)/2 + yn(2:end)/2, ...
         yn(1:end-1).^ry(1)/2 + yn(2:end).^ry(1)/2 ];
     
%         (wn(1:end-1)/2 + wn(2:end)/2).^ry(1)];
 
     
loss  = (eye(n-1) - Theta*pinv(Theta) )*zn;

assignin('base','a12', Theta\zn );

end

