function flux=numflux(ul,ur,param)
s=max(abs(param.fluxd(ul)),abs(param.fluxd(ur)));
flux=0.5*(param.flux(ul)+param.flux(ur)-s*(ur-ul));
end