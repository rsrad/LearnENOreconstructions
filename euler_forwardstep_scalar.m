function unew = euler_forwardstep_scalar(u,h,param)
global non_osc_method
uold=get_boundary_scalar(u,param);
flux=zeros(1, param.nx+1);
f=param.flux(uold);
df=param.fluxd(uold);
% Lax-Friedrichs (LF) flux splitting
a = max(abs(df));
fp = 0.5*(f + a.*uold ); 
fm = 0.5*(f - a.*uold );
switch non_osc_method
    case{3}
        for k=param.gc+1:param.nx+param.gc+1
            fpr = param.reconstruction(fp(k-2),fp(k-1),fp(k),param);
            fmr = param.reconstruction(fm(k+1),fm(k),fm(k-1),param);
            %flux(k-param.gc)=numflux(ul,ur,param);
            flux(k-param.gc)=fpr+fmr;
        end    
    otherwise
        for k=param.gc+1:param.nx+param.gc+1
            fpr = param.reconstruction(fp(k-3),fp(k-2),fp(k-1),fp(k),fp(k+1),param);
            fmr = param.reconstruction(fm(k+2),fm(k+1),fm(k),fm(k-1),fm(k-2),param);
            flux(k-param.gc)=fpr+fmr;
        end
end
unew=u-(h/param.dx)*diff(flux);
end