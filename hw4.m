clear all;
clf;
close all;
randn('seed',1);

beta = zeros(10,1);
beta(3) = 1;
beta(5) = 7;
beta(10) = 3;
global newton_vals;
newton_vals = [];
n=100;
p=10;
mu=10;

X=randn(n,p);
y = X*beta + 0.1*randn(n,1);
lambda = 0.2;

beta = ones(10,1);
t=20*ones(p,1);
x=[beta;t];

delta0=1/lambda;

initial_val=eval_obj(X,y,lambda,p,delta0,x)

global obj_val;
global obj_it;

obj_val = [initial_val];
obj_it = [1];

[opt_x,opt_value]=barrier_lty(X,y,lambda,p,delta0,x,1e-6,mu);

opt_x = opt_x(1:p)

subplot(221);
plot(obj_it,obj_val);
set(gca,'XMinorTick','on','YMinorTick','on');
grid on;
hold on;
plot(obj_it,obj_val,'r*');
xlabel({'Barrier Method Iterations ';'(outer loop)(\mu=10)'}) 
ylabel('Objection Function Value') 
for i=1:size(obj_it,2)
    text(obj_it(i),obj_val(i),['(',(num2str(obj_val(i))),')'],'color','b');
end

subplot(222);
plot(log(newton_vals)/log(10));
grid on;
xlabel({'Newton Method Iterations ';'(inner loop)(\mu=10)'}) 
ylabel('Log10(Duality Gap)') 


%clear all;

newton_vals = [];
mu=100;

delta0=1/lambda;

initial_val=eval_obj(X,y,lambda,p,delta0,x)

obj_val = [initial_val];
obj_it = [1];

[opt_x,opt_value]=barrier_lty(X,y,lambda,p,delta0,x,1e-6,mu);

opt_x = opt_x(1:p)

subplot(223);
plot(obj_it,obj_val);
set(gca,'XMinorTick','on','YMinorTick','on');
grid on;
hold on;
plot(obj_it,obj_val,'r*');
xlabel({'Barrier Method Iterations';' (outer loop)(\mu=100)'}) 
ylabel('Objection Function Value') 
for i=1:size(obj_it,2)
    text(obj_it(i),obj_val(i),['(',(num2str(obj_val(i))),')'],'color','b');
end

subplot(224);
plot(log(newton_vals)/log(10));
grid on;
xlabel({'Newton Method Iterations ';'(inner loop)(\mu=100)'}) 
ylabel('Log10(Duality Gap)') 

set(gcf, 'position', [300 100 1200 800]);