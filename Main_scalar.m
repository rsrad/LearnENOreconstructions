% This matlab code uses trained models and the underlying reconstruction 
% function to reconstruct the ENO3, ENOL3 and WENOC3 states and compute the 
% solution of scalar conservation laws. 
% close all ;
% clear all;
% Define some parameters
global non_osc_method

modelchoice =[1,2,3,4];
% 1=Ms = model trained on smooth data set
% 2= Md = Model trained on discontinuos data set
% 3 = Mu = union model
% 4 = ML = flux limiter

param.test=5; % test case
param.nx=100; % Number of grids
param.rk_method='SSP3'; % ODE solver

param.eps0 = 1.0e-40;            % To use in weno5
param.animation=true;

% % % % % % choice_of_r =1; % 
hidden_layer_n =7;

% non_osc_method variable takes value 1,2,3 for 
% 1: classical ENO3, 2: length based ENO3 and 3: WENOC3 method.
non_osc_method= 3;
NN = true;

% In line function param.reconstruction define the choice of method for
% reconstruction. For reconstruction using model trained on reconsction 
% data set we set it as NN_eno3, NN_enol3 and nn_wenoc3 whereas reconstruction 
% using corresponding subroutines are invoved without prefix "NN_"

if non_osc_method ==1
    method= 'eno3';
    model_dir = 'ENO3/model/';
    prename ='NNENO3-';
    if NN==true
        param.reconstruction= @ NN_eno3;
    else
        param.reconstruction= @ eno3;
    end
    param.gc=4;% no. of ghost cell
    figno=5;
elseif non_osc_method ==2
    method= 'enol3';
    model_dir = 'ENOL3/model/';
    prename ='NNENOL3-';
    param.gc=4;% no. of ghost cell
    if NN==true
        param.reconstruction= @ NN_enol3;
    else
        param.reconstruction= @ enol3;
    end
    figno=6;
else
    method= 'weno3c';
    model_dir = 'WENO3/model/';
    prename ='NNWENO-';
    param.gc=4;% no. of ghost cell
    if NN ==true
        param.reconstruction= @ NN_weno;
    else
        param.reconstruction= @ weno3c;
    end
    figno=7;
end

smooth_weight_file = model_dir +"trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "smoothdata"+".mat";
smooth_bias_file = model_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "smoothdata" + ".mat";

nonsmooth_weight_file = model_dir +"trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "nonsmoothdata"+".mat";
nonsmooth_bias_file = model_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "nonsmoothdata"+ ".mat";

uniondata_weight_file = model_dir +"trained_weights_"+"Hn_"+num2str(hidden_layer_n) + "union_data"+".mat";
uniondata_bias_file = model_dir +"trained_biases_"+"Hn_"+num2str(hidden_layer_n) + "union_data"+ ".mat";


% smooth weights and biases
load(smooth_weight_file,'WEIGHT_smoothdata');
load(smooth_bias_file,'BIAS_smoothdata');
param.weight=WEIGHT_smoothdata;
param.bias=BIAS_smoothdata;

% non_smooth weights and biases
load(nonsmooth_weight_file,'WEIGHT_nonsmoothdata');
load(nonsmooth_bias_file,'BIAS_nonsmoothdata');
param.weight1=WEIGHT_nonsmoothdata; % weight1 and bias1 for non smooth
param.bias1=BIAS_nonsmoothdata;

% union(smooth+discontnuous) weights and biases
load(uniondata_weight_file,'WEIGHT_union_data');
load(uniondata_bias_file,'BIAS_union_data');
param.weight2= WEIGHT_union_data; % weight1 and bias1 for smooth
param.bias2= BIAS_union_data;

for md =1:4
    param.model = modelchoice(md);
    mdl = ' ';
    switch param.model
        case 1
            mdl= 'MS'; % To use Model trained on smooth data set
            mktype ='-or';
            mksize =4;
            intval=1;
            lw=0.5;
            mdd=1;
        case 2
            mdl= 'MD'; % To use Model trained on Discontinuous data set
            mktype ='*b';
            mksize =6;
            intval=1;
            lw=1.0;
            mdd=1;
        case 3
            mdl= 'MU'; % To use Model trained on union data set
            mktype ='sb';
            mksize =6;
            intval=1;
            lw=0.5;
            mdd=3;
        case 4
            mdl= 'ML'; % To use hybrid limiting Model
            mktype ='-ob';
            mksize =6;
            intval=1;
            lw =0.5;   
            mdd=4;
    end
    fntsize=10;
    tic
    [t,u,param] = scalar1D_solver(param); % Main solver
    toc
    %% plot and save
    fprintf('Result is ready. Ploting...\n');
    %%% scheme is in finite volume form
    if  NN==true 
        SchemeName = strcat(prename , mdl , '-Hn', num2str(hidden_layer_n));
       
    else
         str= strsplit(func2str(param.reconstruction), '_');
        SchemeName= string(str(2)); 
    end
    
    
    %'smooth_data','discontinuous_data','Union_data','flux-limiter',NNENOL-ϕ
    %['FV-WENO (',upper(func2str(param.reconstruction)),')'];
    close(2);
    figure(figno);
    plot(param.x,u, mktype,"MarkerSize", mksize,...
        'MarkerIndices',1:intval:length(param.x), "LineWidth", lw,...
        'DisplayName',SchemeName);
    hold on
end
plot(param.x,param.exact_sol,'-k',"LineWidth", 1, 'DisplayName','Exact');
% plot(xex,uex,'-k',"LineWidth", 1, 'DisplayName','Exact');
xlabel('x'); ylabel('u');
axis([-1 1 -0.1 1.2])
%plot(param.x,param.u0,'-r','DisplayName',SchemeName);    
legend('-DynamicLegend','Interpreter', 'none',...
    'Fontsize', fntsize, 'NumColumns',1, 'Location', 'northwest');

