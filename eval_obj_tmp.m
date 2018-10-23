function rs=eval_obj_tmp(X,y,lambda,p,delta,x)

    % object function
    beta = x(1:p);
    t = x(p+1:2*p);

    P=(y-X*beta)'*(y-X*beta);
    Q=lambda*ones(1,p)*abs(beta);

    rs=P+Q;

end