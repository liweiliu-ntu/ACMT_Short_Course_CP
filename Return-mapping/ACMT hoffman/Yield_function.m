function [fnext] = Yield_function(Qnext,Y,b,r)
    fnext=0.5.*(Qnext'*Y*Qnext)+(b'*Qnext).*r-0.5.*r.^2;
end
