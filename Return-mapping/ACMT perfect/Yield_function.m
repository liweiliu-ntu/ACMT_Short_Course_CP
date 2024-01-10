function [fnext] = Yield_function(Qnext,Qy)
    fnext=sqrt(Qnext'*Qnext)-Qy;
end
