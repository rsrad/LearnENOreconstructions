% function for generating classification data for ENO3-length

function [u,r]=weno3c(um1,u0,up1,param)
% global dx
dx =param.dx;
termL = abs(u0 - um1)^2;
termR =abs(u0-up1)^2;
term1= 0.5*(up1-um1)^2;
if ((termL+termR)-term1)^2 < 4*(dx)^4
    w0=1/3;
    r=0;
elseif termL<termR
    w0=1;
    r=1;
else
    w0=0;
    r=2;
end
%if 0.5*(u0+up1)>=0
    u = w0*(1.5*u0 - 0.5*um1) + (1-w0)*0.5*(u0+up1);
% else
%     u = (1-w0)*(1.5*u0 - 0.5*um1) + (w0)*0.5*(u0+up1);


end