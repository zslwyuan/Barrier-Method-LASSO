function rs=eval_obj(X,y,lambda,p,delta,x)

    % object function
    beta = x(1:p);
    t = x(p+1:2*p);

    P=(y-X*beta)'*(y-X*beta);
    Q=lambda*ones(1,p)*t;

    R=0;
    for i=1:p
        R=R+log(t(i)*t(i)-beta(i)*beta(i));
    end
    R=-R/delta;
    rs=P+Q+R;

end