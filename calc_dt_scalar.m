function dt=calc_dt_scalar(u,param)
if param.test==10
    dt=param.dx^(5/3);
else
maxeig=max(abs(param.fluxd(u)));
dt=param.cfl*(param.dx)/maxeig;
end
end