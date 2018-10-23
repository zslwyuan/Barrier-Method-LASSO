function [opt_x,opt_value]=barrier_lty(X,y,lambda,p,delta0,x,error_tol,mu)

    delta = delta0;
    global obj_val;
    global obj_it;
    
    cnt = 1;
    while (1)
        % 1.Centering step.Compute  x'(t) by minimizing tf+φ, subject to Ax=b.
        [newx,newvalue] = backtracking_newton(X,y,lambda,p,delta,x,error_tol);
        
        % 2.Update.x:=x'(t).
        x = newx;  opt_x = x
        opt_value = newvalue    
        obj_val = [obj_val,newvalue];
        cnt = cnt + 1;
        obj_it = [obj_it, cnt];
        
        % 3.Stopping criterion. quit if m/t < error_tol.
        if (2*p/delta<error_tol)
            break
        end
        
        % 4.Increase t.t:=μt
        delta = mu*delta
    end
    
end