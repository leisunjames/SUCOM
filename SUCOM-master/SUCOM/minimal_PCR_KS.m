function [R,t]=minimal_PCR_KS(pts_3d,pts_3d_)

v12=pts_3d(2,:)-pts_3d(1,:);
X_axis=v12'/norm(v12);
v13=pts_3d(3,:)-pts_3d(1,:);
v23=cross(v12,v13);
Y_axis=v23'/norm(v23);
Z_axis=cross(X_axis,Y_axis);

v12=pts_3d_(2,:)-pts_3d_(1,:);
X_axis_=v12'/norm(v12);
v13=pts_3d_(3,:)-pts_3d_(1,:);
v23=cross(v12,v13);
Y_axis_=v23'/norm(v23);
Z_axis_=cross(X_axis_,Y_axis_);

R=[X_axis_,Y_axis_,Z_axis_]*[X_axis,Y_axis,Z_axis]';

p_=sum(pts_3d,1)';
q_=sum(pts_3d_,1)';

p_=p_/3;
q_=q_/3;

t=q_-R*p_;

end