function ul=enol3(um2,um1,u0,up1,up2,param)
global dx
dx=param.dx;
rec_order=3;
c= [1/3,-7/6,11/6;-1/6,5/6,1/3;1/3,5/6,-1/6];
%c= [3/8,-5/4,15/8;-1/8,3/4,3/8;3/8,3/4,-1/8]; %% point value inerpolation  based  ENO reconstruction cofficient(SIAM review shu paper)

state=[um2, um1, u0, up1, up2];

term11 = (u0 - um1)/dx; 
term12 = (2 *u0 - 3 *um1 + um2)/dx;
term13 = um2-2*um1+u0;

term21 = (u0 - up1)/dx; 
term22 = (u0 - um1)/dx;
term23 = um1-2*u0+up1;

term31 = (u0 - up1)/dx; 
term32 = (2 *u0 - 3 *up1 + up2)/dx;
term33 = up2-2*up1+u0;

if (abs(term13)>0)
    beta1= -(1/(2 *(term13)))*dx^2*(term11*sqrt(1 + (term11)^2)...
        - term12* sqrt(1 + (term12)^2) + asinh(term11) - asinh(term12));
else
    beta1=dx *sqrt(1 + term11^2);
end
if (abs(term23)>0)
    beta2= (1/(2 *(term23)))*dx^2*(term21*sqrt(1 + (term21)^2) ...
        + term22* sqrt(1 + (term22)^2) + asinh(term21) + asinh(term22));
else
    beta2=dx *sqrt(1 + term21^2);
end
if (abs(term33)>0)
    beta3= -(1/(2 *(term33)))*dx^2*(term31*sqrt(1 + (term31)^2)...
        - term32* sqrt(1 + (term32)^2) + asinh(term31) - asinh(term32));
else
    beta3=dx *sqrt(1 + term31^2);
end
beta1=beta1.^2;beta2=beta2.^2;beta3=beta3.^2;
L=[beta1;beta2;beta3];
[~,r]=min(L);
r=r-1;
ul=0;
s = rec_order-1;
for i=1:rec_order
    ul=ul+c(r+1,i)*state(r+i); % if r = r-1 labels are 0,1,2
%     ul=ul+c(r,i)*state(r-1+i); % if r =  labels are 1,2,3
%      ul =ul + c(s-r+1,i)*state(s-r+i);            
end

end