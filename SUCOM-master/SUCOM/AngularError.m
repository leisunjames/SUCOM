function Aerr = AngularError(R_opt, R_gt)

Aerr = abs(acos((trace(R_opt'*R_gt)-1)/2));

end