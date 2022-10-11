

%difference for continuous outcome
difference_vec = 0.1:0.05:1;

%difference (in terms of OR) for binary outcomes
difference_vec_bin = linspace(1.1,5,length(difference_vec));

%standard deviation
st_dev = 1;


%logical indicating continuous outcome (0) or binary outcome (1)
model_type = [0, 0, 0, 1, 1, 1];

%trans proportion
trans_prop = 0.03;

%pre-allocation
null_power = zeros(1,length(difference_vec));
power_estimate = zeros(length(difference_vec),length(model_type));
sample_size_vec = zeros(length(difference_vec),length(model_type));

%computing the necessary sample size for each model to maintain a power of
%0.8
for base_model = 1:6
    for ii = 1:length(difference_vec)
        
        %check model type
        if model_type(base_model)
            difference=difference_vec_bin(ii);
        else
            difference=difference_vec(ii);
        end
        
        
        %number of simulations used per iteration of bisection method
        num_runs =100000;

        %run bisection method to compute sample size (per group)
        n_opt = bisection_method(@(x) 0.8-many_hypotheses(x,num_runs,base_model,st_dev,difference,model_type,trans_prop),1,30,1);
        
        %rounding and accounting for the two groups
        n_opt
        sample_size = ceil(n_opt);
        sample_size_vec(ii,base_model) = 2*sample_size; %two because the sample size refers to one group
    end
    save('first_output.mat','difference_vec','difference_vec_bin','sample_size_vec','num_runs')

end

%power calculation under model misspecification
A = zeros(2,2,num_runs);
for test_model = 1:6
    for ii=1:length(difference_vec)
        
        if model_type(test_model)
            difference=difference_vec_bin(ii);
            sample_size = sample_size_vec(ii,4)/2;
        else
            difference=difference_vec(ii);
            sample_size = sample_size_vec(ii,1)/2;
        end

        num_runs =100000;
        if test_model == 1
            %baseline cts outcome
            mu1 = normrnd(0 , st_dev/sqrt(sample_size),[1,num_runs]);
            mu2 = normrnd(difference, st_dev/sqrt(sample_size),[1,num_runs]);

        elseif test_model == 2
            %two-sided misspecification cts
            misspec1 = binornd(1,trans_prop,[ceil(sample_size),num_runs]);
            misspec2 = binornd(1,(1-trans_prop),[ceil(sample_size),num_runs]);

            n1 = normrnd(misspec1*difference,st_dev);
            n2 = normrnd(misspec2*difference,st_dev);

            mu1 = mean(n1,1);
            mu2 = mean(n2,1);

        elseif test_model == 3
            %one-sided misspecification cts
            misspec1 = binornd(1,trans_prop,[ceil(sample_size),num_runs]);

            n1 = normrnd(misspec1*difference,st_dev);

            mu1 = mean(n1,1);
            mu2 = normrnd(difference, st_dev/sqrt(sample_size),[1,num_runs]);
        elseif test_model == 4
            %baseline bin outcome

            %assumes baseline group 1 probability is 0.5
            p_group1=0.5;
            p_group2=((difference*p_group1)/(1-p_group1))/(1+(difference*p_group1)/(1-p_group1));
            
            a=binornd(sample_size, p_group1,[1,1,num_runs]);
            b=binornd(sample_size, p_group2,[1,1,num_runs]);
            
            A(1,1,:) = a;
            A(1,2,:) = b;
            A(2,1,:) = sample_size - a;
            A(2,2,:) = sample_size - b;
        elseif test_model == 5
            %two-sided misspecification bin

            p_group1=0.5;
            p_group2=((difference*p_group1)/(1-p_group1))/(1+(difference*p_group1)/(1-p_group1));
            
            a=binornd(sample_size, p_group1*(1-trans_prop) + p_group2*trans_prop,[1,1,num_runs]);
            b=binornd(sample_size, p_group2*(1-trans_prop) + p_group1*trans_prop,[1,1,num_runs]);
            
            A(1,1,:) = a;
            A(1,2,:) = b;
            A(2,1,:) = sample_size - a;
            A(2,2,:) = sample_size - b;
            
        elseif test_model == 6
            %one-sided misspecification bin

            p_group1=0.5;
            p_group2=((difference*p_group1)/(1-p_group1))/(1+(difference*p_group1)/(1-p_group1));
            
            a=binornd(sample_size, p_group1*(1-trans_prop) + p_group2*trans_prop,[1,1,num_runs]);
            b=binornd(sample_size, p_group2,[1,1,num_runs]);
            

            A(1,1,:) = a;
            A(1,2,:) = b;
            A(2,1,:) = sample_size - a;
            A(2,2,:) = sample_size - b;
            
        end
        
        
        % this is the acceptance formula
        if model_type(test_model)
            accept = zeros(1,num_runs);
            %for kk = 1:num_runs
            parfor kk = 1:num_runs
                accept(kk)=1-fishertest(A(:,:,kk));
            end
        else
            accept=abs(mu1-mu2)<=norminv((1-0.025),0,1)*sqrt((2*st_dev^2)/(sample_size));
        end
        
        % power estimate under the assumed model (just as a validtion step)
        null_power(ii) = 1-normcdf(norminv((1-trans_prop),0,1)*sqrt((2*st_dev^2)/sample_size),difference,sqrt((2*st_dev^2)/sample_size));
        
        % sample based power estimate (if difference is ~= 0)
        power_estimate(ii,test_model)=1-sum(accept)/num_runs;
        
    end
    save('full_output.mat','power_estimate','difference_vec','difference_vec_bin','sample_size_vec','num_runs')

end

save('full_output.mat','power_estimate','difference_vec','difference_vec_bin','sample_size_vec','num_runs')



%bisection method code
function x = bisection_method(func, x1, x2_init, Xtol)
x2=x2_init;
%upper bound search
f1 = func(x1);
f2 = func(x2);
while sign(f1) == sign(f2) && x2<(x2_init*2^10)
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


%code which computes power estimates under different baseline models
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
    accept = zeros(1,num_runs);
    %for kk=1:num_runs
    parfor kk=1:num_runs
     
        accept(kk)=1-fishertest(squeeze(A(:,:,kk)));
    end 
else
    accept=abs(mu1-mu2)<=norminv((1-0.025),0,1)*sqrt((2*st_dev^2)/(ceil(sample_size)));
end


% sample based power estimate (if difference is ~= 0)
power_estimate=1-sum(accept)/num_runs;
end