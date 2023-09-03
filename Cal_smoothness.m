function [r, termL,termR]= Cal_smoothness(um1,u0,up1,dx)

termL = (u0 - um1); 
termR =(u0-up1);

tol = 2*dx^2;
if termL^2 + termR^2< tol
    r = 1;
elseif abs(termR)<eps
    r = (termL/eps)*signum(termR);
else
    r =termL/termR;
end




% if termL^2+ termR^2 <tol
%     r =1;
% elseif termR^2< tol
%     r = (termL/tol)*signum(termR);    
% else
%     r =termL/termR;
% end

% if termL^2 <tol
%     if termR^2> tol
%         r =termL/termR;
%     else
%         r=1.0;        
%     end
% else
%     if termR^2< tol
%        r = (termL/tol)*signum(termR);
%     else
%         r =1;
%     end

end

