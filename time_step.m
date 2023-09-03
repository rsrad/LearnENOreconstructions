function u=time_step(u0,h,param)
switch param.rk_method % options o1, SSP2,SSP3, SSP45
    case 'SSP3'
        % stage 1
        u = euler_forwardstep_scalar(u0,h,param);
        % stage 2
        ui = euler_forwardstep_scalar(u,h,param);
        u = (1/4)*(3*u0+ui);
        % stage 3
        ui = euler_forwardstep_scalar(u,h,param);
        u = (1/3)*(u0+2*ui);
    case 'o1'  
       u=euler_forwardstep_scalar(u0,h,param);  
end