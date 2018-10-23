function [newx,newf_value]=backtracking_newton(X,y,lambda,p,delta,x,error_tol)

    [g,H]=g_H_comp(X,y,lambda,p,delta,x);
    deltax = -inv(H)*g;
    decrement_2 = g'*inv(H)*g;
    t = 0.001;
    newx = x;
    newf_value = eval_obj(X,y,lambda,p,delta,x);
    
    global newton_vals;

    a = 0.1;
    b = 0.9;
    
    % refer to paper: https://web.stanford.edu/~boyd/papers/pdf/l1_ls.pdf
    % the dual value can be obtained:
    s = min(lambda./(abs(2*X'*(X*x(1:p)-y))));
    v = 2*s*(X*x(1:p)-y);

    % 2.Stopping criterion. quit if λ2/2≤ error_tol.
    while (decrement_2/2>error_tol)
        
        % 3.Line search. Choose step size t by backtracking line search.
        newx = x+t*deltax;
        while (eval_obj(X,y,lambda,p,delta,newx) >= eval_obj(X,y,lambda,p,delta,x) +a*t*(g')*deltax)
            t = b*t;
            newx = x+t*deltax;
        end
        
        % 4a.Update. x:=x+t∆xnt
        x = newx;
        newf_value = eval_obj(X,y,lambda,p,delta,x);
        
        % 4b. Calculate the duality gap
        newton_vals = [newton_vals,eval_obj_tmp(X,y,lambda,p,delta,x)-G(v,y)];        
        
        % 1.Compute the Newton step and decrement.
        %    ∆xnt:=−∇^(2)f(x)^(−1)∇f(x);  λ_2:=∇f(x)^(T)∇^(2)f(x)^(−1)∇f(x).
        [g,H]=g_H_comp(X,y,lambda,p,delta,x);
        deltax = -inv(H)*g;
        decrement_2 = g'*inv(H)*g;
    end
    
end

function rs=G(v,y)
    rs = -0.25*v'*v-v'*y;
end