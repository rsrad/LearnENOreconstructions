function r=forward_pass(input,WEIGHT,BIAS) 
% global non_osc_method
nlayer=length(WEIGHT);
% For hidden layers
for i=1:nlayer-1
    output=input*WEIGHT{i}+BIAS{i};
    input=LeakyReLU(output,0.00); %apply activation
end
% For output layers
output=input*WEIGHT{end}+BIAS{end};
output=softmax(output');% apply last activation
[~,ind] = max(output);
r=ind-1;
% if  non_osc_method==1
%     r=ind-1;
% elseif non_osc_method==2
%     r=ind-1;
% else
%     r = ind-1;
% end
end