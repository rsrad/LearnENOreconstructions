function param=initial_conditions_scalar(param)
switch param.test
    case 1
        param.flux=@(u) u;  % flux of linear advection equation
        param.fluxd=@(u) 1; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        %  jumpp=0;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=2;
        param.tstep=param.tf/20;
        param.cfl=0.50;
        
        %  param.u0=(abs(param.x)<=1/3);
                  param.u0= sin(2*pi*param.x);
                  param.exact_sol=@(x,t) sin(2*pi*(param.x-t));
        
        param.u0=sin(pi*param.x).^4;
        param.exact_sol = @ (x,t) (sin(pi*(param.x-t))).^4;
        
        param.bc='Periodic';
    case 2
        param.flux=@(u)0.5*u.^2; % flux of Burgers equation
        param.fluxd=@(u) u; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=0.6;
        param.tstep=param.tf/20;
        param.cfl=0.25;
        % param.u0=0.5+sin(pi*param.x);
        param.u0=-sin(pi*param.x);
        param.bc='Periodic';
    case 3
        param.flux=@(u) 0.5*u.^2; % flux of Burgers equation
        param.fluxd=@(u) u; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=0.5;
        param.tstep=param.tf/20;
        param.cfl=0.5;
        param.u0=(abs(param.x)<=1/3);
        %       param.u0=1*(abs(param.x)<=1/3)-1*(1/3<abs(param.x));
        param.bc='Periodic';
        
        %                 Laney test 4
        b1=-1/3; b2=-(1/3)+ param.tf;
        bshock=(1/3)+(1/2)*param.tf;
        x= param.xmin:param.dx:param.xmax;
        N=length(x);
        for i=1:N
            if x(i)<b1
                uex(i)=0;
            elseif x(i)>=b1 && x(i)<=b2
                uex(i)=(x(i)-b1)/(b2-b1);
            elseif x(i)>b2 && x(i)<=bshock
                uex(i)=1;
            else
                uex(i)=0;
            end
        end
        param.exact_sol =uex;
        
    case 4
        param.flux=@(u) 0.5*u.^2; % flux of Burgers equation
        param.fluxd=@(u) u; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=0.5;
        param.tstep=param.tf/20;
        param.cfl=0.5;
        param.u0=0.5+sin(pi*param.x);
        param.bc='Periodic';
        
    case 5
        param.flux=@(u) u; % flux of linear advection equation
        param.fluxd=@(u) 1; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=2;
        param.tstep=param.tf/20;
        param.cfl=0.6;
        param.u0=(abs(param.x)<=1/3);
        param.exact_sol=param.u0;
        param.bc='Periodic';
    case 6
        param.flux=@(u) 0.5*u.^2; % flux of Burgers equation
        param.fluxd=@(u) u; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=0.3;
        param.tstep=param.tf/20;
        param.cfl=0.5;
        param.u0=1*(abs(param.x)<=1/3)-1*(1/3<abs(param.x));
        param.bc='Periodic';
    case 7
        param.flux=@(u) u;  % flux of linear advection equation
        param.fluxd=@(u) 1; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        %%% jumpp=0;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=2;
        param.tstep=param.tf/20;
        param.cfl=0.50;
        param.u0=(abs(param.x)<=1/3);
        param.bc='Periodic';
        
        
        
    case 8
        param.flux=@(u) u.^2/(u.^2+(1-u).^2); % flux of Buckley leverett equation
        param.fluxd=@(u) (2*u.*(u.^2+(1-u).^2)-u.^2.*(4*u-2))/(u.^2 +(1-u).^2).^2; % flux derivative
        param.xmin=-1;
        param.xmax=1;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=2.0;
        param.tstep=param.tf/20;
        param.cfl=0.5;
        param.u0=-1*(param.x<=0)+0*(param.x>0);
        param.bc='Periodic';
        
        
    case 9
        
        param.flux=@(u) 0.5*u.^2; % flux of Burgers equation
        param.fluxd=@(u) u; % flux derivative
        param.xmin=0;
        param.xmax=1.4;
        param.dx=(param.xmax-param.xmin)/(param.nx-1);
        param.x=param.xmin:param.dx:param.xmax;
        param.tf=1.4;
        param.tstep=param.tf/20;
        param.cfl= 0.5;
        param.u0= 10.*(param.x-0.2).*(0.2<param.x<=0.3)+ 10.*(0.4-param.x).*(0.3<param.x<=0.4)+1*(0.6<param.x<=0.8)+100.*(param.x-1).*(1.2-param.x).*(0.2<param.x<=0.3)+0*(param.x>0);
        param.bc='Periodic';
        
        
end
end