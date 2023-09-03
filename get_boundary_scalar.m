function ubc=get_boundary_scalar(u,param)
nx=numel(u);
ubc=zeros(1,nx+2*param.gc);
switch param.bc
    case 'Neumann'
            ubc(1:param.gc)=u(1);
            ubc(param.gc+1:nx+param.gc)=u;
            ubc(nx+param.gc+1:nx+2*param.gc)=u(end);
     case 'Periodic'
            ubc(1:param.gc)=u(nx-param.gc:nx-1);
            ubc(param.gc+1:param.gc+nx)=u;
            ubc(nx+param.gc+1:nx+2*param.gc)=u(2:param.gc+1);
        %end
end