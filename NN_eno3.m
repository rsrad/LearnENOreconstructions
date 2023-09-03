function ul=NN_eno3(um2,um1,u,up1,up2,param)
% NN_eno(um2,um1,u,up1,up2,dx) and use global param to run reconstruction
% accuracy mat file code.
c= [1/3,-7/6,11/6;-1/6,5/6,1/3;1/3,5/6,-1/6];
dx = param.dx;
input=[um2, um1, u, up1, up2];

rec_order=3;
s = rec_order-1;
ul_smooth=0; ul_jump=0; ul_union=0;
m = mean(input); sd=std(input); 
standardize_input= (input-m)/(sd+10^-35);
standardize_input= [standardize_input,dx];

switch param.model
    case{1}
        %%Evaluating the stencil shift by NN model 1
        r1=forward_pass(standardize_input,param.weight,param.bias); 
         for i=1:rec_order            
            ul_smooth =ul_smooth+ c(s-r1+1,i)*input(s-r1+i);  
        end
        ul = ul_smooth;        
        
    case{2}
        % Evaluating the stencil shift by NN model2 ;
        r2=forward_pass(standardize_input,param.weight1,param.bias1);
        for i=1:rec_order
            ul_jump =ul_jump+ c(s-r2+1,i)*input(s-r2+i);
        end
        ul= ul_jump;
    case{3}
        % Evaluating the stencil shift by NN model3 ;
        r3=forward_pass(standardize_input,param.weight2,param.bias2);
        for i=1:rec_order
            ul_union =ul_union+ c(s-r3+1,i)*input(s-r3+i);
        end
        ul= ul_union;
    case{4}
        %  model limiter based ENO reconstructed states
        r1=forward_pass(standardize_input,param.weight,param.bias);
        r2=forward_pass(standardize_input,param.weight1,param.bias1);
        for i=1:rec_order
            ul_smooth =ul_smooth+ c(s-r1+1,i)*input(s-r1+i);
            ul_jump =ul_jump+ c(s-r2+1,i)*input(s-r2+i);
        end
      % function phi
        [r]= Cal_smoothness(um1,u,up1,param.dx);
        ptr =3.0;
        phi = abs(r)>1/ptr && abs(r)<ptr;
        ul=(phi*ul_smooth+(1-phi)*ul_jump);        
end
end