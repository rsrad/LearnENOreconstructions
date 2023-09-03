function ul=NN_weno(um1,u,up1,param)
% NN_eno(um2,um1,u,up1,up2,dx) and use global param to run reconstruction
% accuracy mat file code.
c= [1/3,-7/6,11/6;-1/6,5/6,1/3;1/3,5/6,-1/6];
%c= [3/8,-5/4,15/8;-1/8,3/4,3/8;3/8,3/4,-1/8]; %% point value inerpolation  based  ENO reconstruction rcofficient(SIAM review shu paper)
% global param
dx = param.dx;
input=[um1, u, up1];
m = mean(input);
sd=std(input);
standardize_input= (input-m)/(sd+10^-35);
standardize_input= [standardize_input,dx];
switch param.model
    case{1}
        % To use trained network use forward pass
        r1=forward_pass(standardize_input,param.weight,param.bias);
        if r1 ==1
            w0=1/3;
        elseif r1==0
            w0=1;
        else
            w0=0;
        end
        ul = w0*(1.5*u - 0.5*um1) + (1-w0)*0.5*(u+up1);       
    case{2}
        r1=forward_pass(standardize_input,param.weight1,param.bias1);
        if r1 ==1
            w0=1/3;
        elseif r1==0
            w0=1;
        else
            w0=0;
        end
        ul = w0*(1.5*u - 0.5*um1) + (1-w0)*0.5*(u+up1);
    case{3}
        r1=forward_pass(standardize_input,param.weight2,param.bias2);
        if r1 ==1
            w0=1/3;
        elseif r1==0
            w0=1;
        else
            w0=0;
        end
        ul = w0*(1.5*u - 0.5*um1) + (1-w0)*0.5*(u+up1);
    case{4}
        r1=forward_pass(standardize_input,param.weight,param.bias);
        if r1 ==1
            w0=1/3;
        elseif r1==0
            w0=1;
        else
            w0=0;
        end
        ul = w0*(1.5*u - 0.5*um1) + (1-w0)*0.5*(u+up1);
        r2=forward_pass(standardize_input,param.weight1,param.bias1);
        if r2 ==1
            w0=1/3;
        elseif r2==0
            w0=1;            
        else
            w0=0;
        end        
        ul1 = w0*(1.5*u - 0.5*um1) + (1-w0)*0.5*(u+up1);       
%         Switch function phi        
        [r]= Cal_smoothness(um1,u,up1,param.dx);
        ptr =3.0;
        phi = abs(r)>1/ptr && abs(r)<ptr;       
        ul= phi*ul+(1-phi)*ul1;
  
end
end
