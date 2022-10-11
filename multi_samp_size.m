%simlar script to power_calc_script but this gives repeated calculations of
%the sample size (so that error bounds may be computed)

difference_vec = 0.1:0.05:1;
difference_vec_bin = linspace(1.1,5,length(difference_vec));
st_dev = 1;
base_model = 1;
test_model = 2;

model_type = [0, 0, 0, 1, 1, 1];
trans_prop = 0.03;

num_samps = 30;
sample_size_vec = zeros(length(difference_vec),length(model_type),num_samps);

options = optimset('TolX',1);

for base_model = 1:6
    for ii = 1:length(difference_vec)
        for sample=1:num_samps
             if model_type(base_model)
                 difference=difference_vec_bin(ii);
             else
                 difference=difference_vec(ii);
             end
      
            num_runs =10000;
            
            x1 = 1;
            x2 = 30;
            n_opt = bisection_method(@(x) 0.8-many_hypotheses(x,num_runs,base_model,st_dev,difference,model_type,trans_prop),x1,x2,1);
            
            n_opt
            sample_size = ceil(n_opt);
            sample_size_vec(ii,base_model,sample) = 2*sample_size; %two because the sample size refers to one group
        end
        save('multisamp_output1.mat','difference_vec','difference_vec_bin','sample_size_vec','num_runs')
    end
end

%computing the mean and variance of the proprotional difference in sample size estimate compared to baseline
mean_base1=mean(100*((sample_size_vec(:,[2,3,5,6],:)-sample_size_vec(:,[1,1,4,4],:))./sample_size_vec(:,[1,1,4,4],:)),3);
var_base1 = var(100*((sample_size_vec(:,[2,3,5,6],:)-sample_size_vec(:,[1,1,4,4],:))./sample_size_vec(:,[1,1,4,4],:)),0,3);

%error bounds
l_b = mean_base1-tinv(0.975,num_samps-1).*(sqrt(var_base1/num_samps));
u_b = mean_base1+tinv(0.975,num_samps-1).*(sqrt(var_base1/num_samps));


%computing the mean sample size estimates and error bounds for the baseline
%models
mean_cts = mean(sample_size_vec(:,1,:),3);
mean_bin = mean(sample_size_vec(:,4,:),3);


log_m_cts = log(mean_cts);
log_m_bin = log(mean_bin);

ub_cts = log(mean_cts + tinv(0.975,num_samps-1).*(sqrt(var(sample_size_vec(:,1,:),0,3)/num_samps)));
ub_bin = log(mean_bin + tinv(0.975,num_samps-1).*(sqrt(var(sample_size_vec(:,2,:),0,3)/num_samps)));


lb_cts = log(mean_cts - tinv(0.975,num_samps-1).*(sqrt(var(sample_size_vec(:,1,:),0,3)/num_samps)));
lb_bin = log(mean_bin - tinv(0.975,num_samps-1).*(sqrt(var(sample_size_vec(:,2,:),0,3)/num_samps)));



function x = bisection_method(func, x1, x2_init, Xtol)
x2=x2_init;
%upper bound search
f1 = func(x1);
f2 = func(x2);
while sign(f1) == sign(f2) && x2<min((x2_init*2^10),100000)
    f1 = f2;
    x2 = 2*x2;
    f2 = func(x2);
    if x2>=(x2_init*2^10)
        display('big samp')
    end
end

%bisection
while abs(x1-x2)>Xtol

    f_new = func((x1+x2)/2);

    if sign(f_new)~=sign(f2)
        x1=(x1+x2)/2;  
    else
        x2 = (x1 + x2)/2;
        f2 = f_new;

    end
end
x = x2;
end


function [power_estimate]=many_hypotheses(sample_size,num_runs,model,st_dev,difference,model_type,trans_prop)

if model ==1
    mu1 = normrnd(0 , st_dev/sqrt(ceil(sample_size)),[1,num_runs]);
    mu2 = normrnd(difference, st_dev/sqrt(ceil(sample_size)),[1,num_runs]);
elseif model ==2

    misspec1 = binornd(1,trans_prop,[ceil(sample_size),num_runs]);
    misspec2 = binornd(1,(1-trans_prop),[ceil(sample_size),num_runs]);

    n1 = normrnd(misspec1*difference,st_dev);
    n2 = normrnd(misspec2*difference,st_dev);

    mu1 = mean(n1,1);
    mu2 = mean(n2,1);
elseif model == 3

    misspec1 = binornd(1,trans_prop,[ceil(sample_size),num_runs]);

    n1 = normrnd(misspec1*difference,st_dev);
    
    mu1 = mean(n1,1);
    mu2 = normrnd(difference, st_dev/sqrt(sample_size),[1,num_runs]);

elseif model == 4
    
    p_group1=0.5;
    p_group2=((difference*p_group1)/(1-p_group1))/(1+(difference*p_group1)/(1-p_group1));
    
    a=binornd(ceil(sample_size), p_group1,[1,1,num_runs]);
    b=binornd(ceil(sample_size), p_group2,[1,1,num_runs]);
    
    A = zeros(2,2,num_runs);
    
    A(1,1,:) = a;
    A(1,2,:) = b;
    A(2,1,:) = ceil(sample_size) - a;
    A(2,2,:) = ceil(sample_size) - b;
    
elseif model == 5
    
    
    p_group1=0.5;
    p_group2=((difference*p_group1)/(1-p_group1))/(1+(difference*p_group1)/(1-p_group1));

    a=binornd(ceil(sample_size), p_group1*(1-trans_prop) + p_group2*trans_prop,[1,1,num_runs]);
    b=binornd(ceil(sample_size), p_group2*(1-trans_prop) + p_group1*trans_prop,[1,1,num_runs]);

    A = zeros(2,2,num_runs);
    
    A(1,1,:) = a;
    A(1,2,:) = b;
    A(2,1,:) = ceil(sample_size) - a;
    A(2,2,:) = ceil(sample_size) - b;
    
elseif model == 6

    p_group1=0.5;
    p_group2=((difference*p_group1)/(1-p_group1))/(1+(difference*p_group1)/(1-p_group1));

    a=binornd(ceil(sample_size), p_group1*(1-trans_prop) + p_group2*trans_prop,[1,1,num_runs]);
    b=binornd(ceil(sample_size), p_group2,[1,1,num_runs]);
    

    A = zeros(2,2,num_runs);
    
    A(1,1,:) = a;
    A(1,2,:) = b;
    A(2,1,:) = ceil(sample_size) - a;
    A(2,2,:) = ceil(sample_size) - b;
    

end

% this is the acceptance formula
if model_type(model)
    display('hello joe')
    accept = zeros(1,num_runs);
    %for kk=1:num_runs
    parfor kk=1:num_runs
     
        accept(kk)=1-fishertest(squeeze(A(:,:,kk)));
    end 
else
    accept=abs(mu1-mu2)<=norminv((1-0.025),0,1)*sqrt((2*st_dev^2)/(ceil(sample_size)));
end

power_estimate=1-sum(accept)/num_runs;
end
