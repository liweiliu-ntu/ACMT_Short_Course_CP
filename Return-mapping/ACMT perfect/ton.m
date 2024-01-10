function [t_on,Qt,qt,d_ton] = ton(Q_dot1,Q_off,q_dot1,tton, q_off,Qy)
syms x
eqn=sqrt((Q_off+x.*Q_dot1)'*(Q_off+x.*Q_dot1))-Qy==0;
        d_ton=double(max(solve(eqn,x)));
        t_on=double(max(solve(eqn,x)))+tton;
         Qt=Q_off+Q_dot1.*(t_on-tton);
         qt=q_off+(q_dot1)*d_ton;
         
end