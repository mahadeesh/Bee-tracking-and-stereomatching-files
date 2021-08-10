function [XL_MID_PT,XR_MID_PT,X2C] = stereo_triangulation_camera_positions(xL,xR,om,T,fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,kc_right,alpha_c_right)

% This function computes left and right cameras 3D position in the world coordinates. 

% This function was adapted from (c) Jean-Yves Bouguet - Intel Corporation - April 9th, 2003, by mandiyam mahadeeswara(m.mandiyam@uq.edu.au)
%
% [XL,XR] = stereo_triangulation(xL,xR,om,T,fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,kc_right,alpha_c_right),
%
% Function that computes the position of a set on N points given the left and right image projections.
% The cameras are assumed to be calibrated, intrinsically, and extrinsically.
%
% Input:
%           xL: 2xN matrix of pixel coordinates in the left image
%           xR: 2xN matrix of pixel coordinates in the right image
%           om,T: rotation vector and translation vector between right and left cameras (output of stereo calibration)
%           fc_left,cc_left,...: intrinsic parameters of the left camera  (output of stereo calibration)
%           fc_right,cc_right,...: intrinsic parameters of the right camera (output of stereo calibration)
%
% Output:
%
%           XL: 3xN matrix of coordinates of the points in the left camera reference frame
%           XR: 3xN matrix of coordinates of the points in the right camera reference frame
%
% Note: XR and XL are related to each other through the rigid motion equation: XR = R * XL + T, where R = rodrigues(om)
% For more information, visit http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
%


%--- Normalize the image projection according to the intrinsic parameters of the left and right cameras
xt = norm_img_proj(xL,fc_left,cc_left,kc_left,alpha_c_left);
xtt = norm_img_proj(xR,fc_right,cc_right,kc_right,alpha_c_right);

%--- Extend the normalized projections in homogeneous coordinates
xt = [xt;ones(1,size(xt,2))];
xtt = [xtt;ones(1,size(xtt,2))];

%--- Number of points:
N = size(xt,2);

%--- Rotation matrix corresponding to the rigid motion between left and right cameras:
R = rodrigues(om);


%--- Triangulation of the rays in 3D space:

u = R * xt;

n_xt2 = dot(xt,xt); % u.u
n_xtt2 = dot(xtt,xtt); % v.v 

T_vect = repmat(T, [1 N]);

DD = n_xt2 .* n_xtt2 - dot(u,xtt).^2;  % ac - b^2

dot_uT = dot(u,T_vect);  % u.w0 
dot_xttT = dot(xtt,T_vect); % v.w0
dot_xttu = dot(u,xtt); % u.v

NN1 = dot_xttu.*dot_xttT - n_xtt2 .* dot_uT;
NN2 = n_xt2.*dot_xttT - dot_uT.*dot_xttu;

Zt = NN1./DD;
Ztt = NN2./DD;

X1 = xt .* repmat(Zt,[3 1]);
X2 = R'*(xtt.*repmat(Ztt,[3,1])  - T_vect);

Err = norm(X1) - norm(X2);

X2C = R' * ( -T_vect);

XL_MID_PT = X1;

XR_MID_PT = R'*(xtt.*repmat(Ztt,[3,1])) ;





