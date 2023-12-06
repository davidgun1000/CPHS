function [B_var] = jitChol(B_var)
%keyboard
[R,p]=chol(B_var);
if p>0
    min_eig=min(eig(B_var));
    d=size(B_var,1);
    delta=max(0,-2*min_eig+10^(-5)).*eye(d);
    B_var=B_var+delta;
else
    B_var=B_var;
end

end
