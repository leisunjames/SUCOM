function [s_raw,R_raw,t_raw]=minimal_PCR_UKS(pts_3d,pts_3d_,noise)

v12=pts_3d((2),:)-pts_3d((1),:);
X_axis=v12'/norm(v12);
v13=pts_3d((3),:)-pts_3d((1),:);
v23=cross(v12,v13);
Y_axis=v23'/norm(v23);
Z_axis=cross(X_axis,Y_axis);

v12=pts_3d_((2),:)-pts_3d_((1),:);
X_axis_=v12'/norm(v12);
v13=pts_3d_((3),:)-pts_3d_((1),:);
v23=cross(v12,v13);
Y_axis_=v23'/norm(v23);
Z_axis_=cross(X_axis_,Y_axis_);

R_raw=[X_axis_,Y_axis_,Z_axis_]*[X_axis,Y_axis,Z_axis]';

p_=sum(pts_3d,1)';
q_=sum(pts_3d_,1)';

p_=p_/3;
q_=q_/3;

s_=zeros(3,1);s_scale=zeros(3,1);s__=zeros(3,1);

for i=1:3
    s_(i)=norm(pts_3d_((i),:)'-q_)/norm(pts_3d((i),:)'-p_);
    s_scale(i)=1/((7*noise)^2/(norm(pts_3d((i),:)'-p_))^2);
    s__(i)=s_(i)*s_scale(i);
end

s_raw=sum(s__)/sum(s_scale);

t_raw=q_-s_raw*R_raw*p_;

end