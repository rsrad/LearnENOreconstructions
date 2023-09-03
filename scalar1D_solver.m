function [t,u,param]=scalar1D_solver(param)
param=initial_conditions_scalar(param);
% fprintf('Solver is starting with the following parameters:\n')
% disp(param)
t=0;
u=param.u0;
while t<param.tf
    dt=calc_dt_scalar(u,param);
    if (t+dt > param.tf)
        dt = param.tf - t;
    end
    u=time_step(u,dt,param);
    t=t+dt;
    if param.animation
        figure(2)
        plot(u,'-s')
        drawnow
    end
end
end