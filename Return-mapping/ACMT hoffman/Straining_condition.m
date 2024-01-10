function [Sg] = Straining_condition(Y,K,r,b,q_dot1,Qnow)


Sg=Qnow'*Y*K*q_dot1+r*b'*K*q_dot1;
end