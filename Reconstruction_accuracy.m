clear all
close all 
% Define some parameters
global param
%BC='Periodic';
order='3rd_Reconst';
fun_choice =1;  cell_average_function = 1;

    %%%%%%%
modelchoice =[1,2,3,4,5];
param.model = modelchoice(1);
choice_of_r =1;
hidden_layer_n =7;

% 1= Ms = Model trained on smooth data set
% 2= Md = Model trained on discontinuos data set
% 3 = Mu = union model
% 4 = Mh = hybrid model
% 5 = ML = flux limiter

eno_method= 2; %1: classical, 2: length based eno
if eno_method == 1
    method= 'eno3';
    data_dir = '../ENO3data/Details/';
else
    method= 'enol3';
    data_dir = '../ENOL3data/Details/';
end

smooth_weight_file = data_dir + "trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "smoothdata_"+choice_of_r+".mat";
smooth_bias_file = data_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "smoothdata_" + choice_of_r+ ".mat";

nonsmooth_weight_file = data_dir +"trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "nonsmoothdata_"+choice_of_r+".mat";
nonsmooth_bias_file = data_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "nonsmoothdata_" + choice_of_r+ ".mat";

uniondata_weight_file = data_dir +"trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "union_data_"+choice_of_r+".mat";
uniondata_bias_file = data_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "union_data_" + choice_of_r+ ".mat";

% smooth weights and biases
 load(smooth_weight_file,'WEIGHT_smoothdata')
 load(smooth_bias_file,'BIAS_smoothdata')
 param.weight=WEIGHT_smoothdata;
 param.bias=BIAS_smoothdata;

% non_smooth weights and biases
 load(nonsmooth_weight_file,'WEIGHT_nonsmoothdata')
 load(nonsmooth_bias_file,'BIAS_nonsmoothdata')
 param.weight1=WEIGHT_nonsmoothdata; % weight1 and bias1 for non smooth
 param.bias1=BIAS_nonsmoothdata; 

% union(smooth+discontnuous) weights and biases
 load(uniondata_weight_file,'WEIGHT_union_data')
 load(uniondata_bias_file,'BIAS_union_data')
 param.weight2= WEIGHT_union_data; % weight1 and bias1 for smooth
 param.bias2= BIAS_union_data;

%%%%%%%%%%%%%
switch fun_choice
    case 1
        fun =@(x1)sin(pi*x1);
        xlb=-1;xrb=1;
    case 2
        fun=@(x1)(x1<=0);
        xlb=-1;xrb=1;
    case 3
        fun = @(x1)sin(pi*x1).^4;
        xlb=-1;xrb=1;
    case 4
        fun=@(x1)(x1<=0);
        xlb=-1;xrb=1;
end
export_result_as = 'PrintedLatex';

%%
if strcmp(export_result_as,'textfile')
     %    delete 'err.txt'
    fid = fopen('err.txt', 'w+');
    fprintf(fid, ' N 	      L^{\\infty} error 	 Rate     	 L^1 error     	 Rate \n');
    %     fclose(fid);
end
if strcmp(export_result_as,'textPrinted')
    fprintf( ' N        L^{\\infty} error   Rate       L^1 error 	   	 Rate \n');
end
for n=2:8
    % Number of grid points in domain [xlb, xrb]
    N=10*2^(n-1)-1;
    %% Main Code
    x=linspace(xlb,xrb,N); % Equal space x_{1},x_{2},..., x_{N}
    dx=x(3)-x(2);
    dx1= ones(1,N)*dx;
    ext_x = [xlb-2*dx, xlb-dx, x, xrb+dx, xrb+2*dx];
    ext_N=length(ext_x);
   % u=fun(ext_x);
    
    [il2, il1, ii,ir1,ir2]=apply_bc_indices(ext_N);
    xhalf_p=(ext_x(ii)+ ext_x(ir1))/2;
    xhalf_m= xhalf_p-dx;
    x3half_p =xhalf_p+dx;
    x3half_m =xhalf_m-dx;
    x5half_p= x3half_p-dx;
    x5half_m= x3half_m-dx;
    switch cell_average_function
        case 1
            %%%%%   cell average value based ENO reconstruction for sin(pi*x)

            um2 = ((1/dx).*(-1/pi).*(cos(pi*x3half_m)-cos(pi*x5half_m)));
            um1 = (1/dx).*((-1/pi).*(cos(pi*xhalf_m)-cos(pi*x3half_m)));
            u0 = (1/dx).*((-1/pi).*(cos(pi*xhalf_p)-cos(pi*xhalf_m)));
            up1 = (1/dx).*((-1/pi).*(cos(pi*x3half_p)-cos(pi*xhalf_p)));
            up2 = (1/dx).*((-1/pi).*(cos(pi*x5half_p)-cos(pi*x3half_p)));
        case 2
            %%%%% cell average value based ENO reconstruction for sin(2*pi*x)
            um2 = (1/dx).*((-1/(2*pi)).*(cos(2*pi*x3half_m)-cos(2*pi*x5half_m)));
            um1 = (1/dx).*((-1/(2*pi)).*(cos(2*pi*xhalf_m)-cos(2*pi*x3half_m)));
            u0 = (1/dx).*((-1/(2*pi)).*(cos(2*pi*xhalf_p)-cos(2*pi*xhalf_m)));
            up1 = (1/dx).*((-1/(2*pi)).*(cos(2*pi*x3half_p)-cos(2*pi*xhalf_p)));
            up2 = (1/dx).*((-1/(2*pi)).*(cos(2*pi*x5half_p)-cos(2*pi*x3half_p)));
        case 3
            %%%%%%%% cell average value based ENO reconstruction for sin(pix)^4
            um2= (1/dx).*((3*x3half_m)/8 - (1/pi).*(sin(2*pi*x3half_m)/4 - sin(4*pi*x3half_m)/32)-  ((3*x5half_m)/8 - (1/pi).*(sin(2*pi*x5half_m)/4 - sin(4*pi*x5half_m)/32)));
            um1= (1/dx).*((3*xhalf_m)/8 - (1/pi).*(sin(2*pi*xhalf_m)/4 - sin(4*pi*xhalf_m)/32)-  ((3*x3half_m)/8 - (1/pi).*(sin(2*pi*x3half_m)/4 - sin(4*pi*x3half_m)/32)));
            u0= (1/dx).*((3*xhalf_p)/8 - (1/pi).*(sin(2*pi*xhalf_p)/4 - sin(4*pi*xhalf_p)/32)-  ((3*xhalf_m)/8 - (1/pi).*(sin(2*pi*xhalf_m)/4 - sin(4*pi*xhalf_m)/32)));
            up1= (1/dx).*((3*x3half_p)/8 - (1/pi).*(sin(2*pi*x3half_p)/4 - sin(4*pi*x3half_p)/32)-  ((3*xhalf_p)/8 - (1/pi).*(sin(2*pi*xhalf_p)/4 - sin(4*pi*xhalf_p)/32)));
            up2= (1/dx).*((3*x5half_p)/8 - (1/pi).*(sin(2*pi*x5half_p)/4 - sin(4*pi*x5half_p)/32)-  ((3*x3half_p)/8 - (1/pi).*(sin(2*pi*x3half_p)/4 - sin(4*pi*x3half_p)/32)));
    end
    %%%%% grid value based ENO reconstruction

%     um2 = fun(ext_x(il2));
%     um1 = fun(ext_x(il1));
%     u0 = fun(ext_x(ii));
%     up1 = fun(ext_x(ir1));
%     up2 = fun(ext_x(ir2));

%size(param.weight{1,1})
%size(um2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     um2 = (1/dx).*((-1/pi).*(cos(pi*x3half_m)-cos(pi*x5half_m)));
%     um1 = (1/dx).*((-1/pi).*(cos(pi*xhalf_m)-cos(pi*x3half_m)));
%     u0 = (1/dx).*((-1/pi).*(cos(pi*xhalf_p)-cos(pi*xhalf_m)));
%     up1 = (1/dx).*((-1/pi).*(cos(pi*x3half_p)-cos(pi*xhalf_p)));
%     up2 = (1/dx).*((-1/pi).*(cos(pi*x5half_p)-cos(pi*x3half_p)));
    %% %%%%%%%%%%%%%%%%%%%%%%%%
%size(param.weight{1,1})
%size(um2)
[ul]=arrayfun(@NN_eno,um2,um1,u0,up1,up2,dx1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch order
        case '3rd_Reconst'
            Recons_uhalf=ul;
    end

    s = fun(xhalf_p);
    
    %% *******************************************************
     ein(n)=max(abs(s- Recons_uhalf));
     el1(n)=min(dx)*sum(abs(s- Recons_uhalf));
    
    if n==1
         err_rate_in(n)=0.0;
    else
         err_rate_in(n)=log2(ein(n-1)/ein(n)); %/log(delx(n)/delx(n-1));
    end
    if n==1
         err_rate_l1(n)=0.0;
        if strcmp(export_result_as,'PrintedLatex')
                        delete 'err.txt'
            disp('%% Exporting latex table format of accuracy test')
            disp('\begin{table}[htb!]')
            disp('\centering')
            disp('\begin{tabular}{|c|c|c|c|c|}')
            fprintf('\\hline N  & $L^\\infty$ error  &  Rate &  $L^1$ error & Rate \\\\ \n')
            fprintf('\\hline %d & %.5f &  ...  & %.5f & ... \\\\ \n', N+1, ein(n), el1(n))
           
        end
    else
         err_rate_l1(n)=log2(el1(n-1)/el1(n)); %/log(delx(n)/delx(n-1));
        if strcmp(export_result_as,'PrintedLatex')
            fprintf('\\hline %d & %.5e &  %.2f & %.5e & %.2f \\\\ \n',N+1,ein(n), err_rate_in(n), el1(n), err_rate_l1(n))
        end
    end
    
    if strcmp(export_result_as,'textfile')
        %         fid = fopen('err.txt', 'a+');
        if N+1<100
            fprintf(fid, ' %d \t %.5e \t %.2f \t %.5e \t %.2f \n',N+1,ein(n), err_rate_in(n), el1(n), err_rate_l1(n));
        else
            fprintf(fid, '%d \t %.5e \t %.2f \t %.5e \t %.2f \n',N+1,ein(n), err_rate_in(n), el1(n), err_rate_l1(n));
        end
        %         fclose(fid);
    end
    if strcmp(export_result_as,'textPrinted')
        if N+1<100
            fprintf( ' %d \t %.5e \t %.2f \t %.5e \t %.2f \n',N+1,ein(n), err_rate_in(n), el1(n), err_rate_l1(n));
        else
            fprintf( '%d \t %.5e \t %.2f \t %.5e \t %.2f \n',N+1,ein(n), err_rate_in(n), el1(n), err_rate_l1(n));
        end
    end
end
if strcmp(export_result_as,'PrintedLatex')
    disp('\hline ')
    disp('\end{tabular} ')
    disp('\caption{}')
    disp('\label{tab:}')
%    fprintf('%% The solution of test %d at T=%.2f for CFL=%.2f. \n',u,t,CFL)
    disp(['% It was exucuted through ',cd,' on ',date])
    disp('\end{table}')
end