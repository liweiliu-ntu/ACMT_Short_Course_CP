function [Sg] = Straining_condition(q_dot1,Qnow)
Sg=Qnow'*q_dot1;
end