close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k)';
 % this will be useful in a bit
%%% set up the discrete value of A
a1 = 1.1;
a2 = 0.678;
A= [a1;a2];
num_a = 2;
%%% set up the transition matrix
ph = 0.977;
p_h = 1 - ph;
pl = 0.926;
p_l = 1-pl;
p = [ph p_h;p_l pl];

%%% _
cons = zeros(num_k,num_a);
ret = zeros (num_k,num_a);
v_guess = zeros(num_k,num_a);
value_mat = zeros(num_k,num_a);
vfn = zeros(num_k,num_a);
s = zeros (num_k,num_a);
%%%% Iteration
dis = 1; tol = 1e-06;iter = 1; % tolerance for stopping 
while dis > tol
    for i = 1:num_k
        for j = 1:num_a;
            cons = A(j)*k(i)^alpha+(1-delta)*k(i)-k;
            neg = find(cons<0);
            cons(neg) = -Inf;
            ret(:,j)=(cons.^(1-sigma)-1)/(1-sigma);
            ret(neg,j)=-inf;
        end
        [vfn(i,:),dr(i,:)] = max(ret+beta*(v_guess*p));
    end;
       dis = max(abs(vfn - v_guess));
       v_guess = vfn;
       iter = iter + 1 ;
end
kp =  k(dr);
k_mat = [k k];
s= kp-(1-delta)*k_mat;
plot(k,vfn);
figure;
plot(k,kp);
figure;
plot(k,s);


            
