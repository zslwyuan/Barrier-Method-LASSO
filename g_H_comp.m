function [g,H]=g_H_comp(X,y,lambda,p,delta,x)

    beta = x(1:p);
    t = x(p+1:2*p);

    %compute the gradient
    g_b_0 = 2*X'*(X*beta-y);
    g_b_1 = zeros(p,1);
    for i=1:p
        g_b_1(i)=2*beta(i)/(t(i)*t(i)-beta(i)*beta(i));
    end
    g_b_1 = g_b_1/delta;
    g_b = g_b_0 + g_b_1;

    g_t_0 = lambda*ones(p,1);
    g_t_1 = zeros(p,1);
    for i=1:p
        g_t_1(i)=2*t(i)/(t(i)*t(i)-beta(i)*beta(i));
    end
    g_t_1 = -g_t_1/delta;
    g_t = g_t_0 + g_t_1;

    g = [g_b;g_t];

    %compute the hessian matrix
    P1_v = zeros(p,1);
    for i=1:p
        P1_v(i) = 2*(t(i)^2+beta(i)^2)/(t(i)^2-beta(i)^2)^2/delta;
    end
    P1 = diag(P1_v);

    P2_v = zeros(p,1);
    for i=1:p
        P2_v(i) = -4*(t(i)*beta(i))/(t(i)^2-beta(i)^2)^2/delta;
    end
    P2 = diag(P2_v);

    H = [2*X'*X+P1, P2;P2,P1];

end