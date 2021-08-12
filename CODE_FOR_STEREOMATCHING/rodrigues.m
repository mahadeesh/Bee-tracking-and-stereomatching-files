function	[OUT,dOUT]=rodrigues(IN)

% RODRIGUES rotation formula.
%
% 		If IN is a 3x3 rotation matrix then OUT is the
%		corresponding 3x1 rotation vector
% 		if IN is a rotation 3-vector then OUT is the
%		corresponding 3x3 rotation matrix
%

% Copyright (c)     1993    Pietro Perona California Institute of Technology
%                   2003    Jean-Yves Bouguet
%                   2007    Mike Burl

[m,n] = size(IN);

bigeps = 10e+20*eps;

%input is rotational vector
if ((m==1) && (n==3)) || ((m==3) && (n==1))
    theta = norm(IN);
    if theta < eps
        
        R = eye(3);
        dRdIN = [0 0 0;0 0 1;0 -1 0;0 0 -1;0 0 0;1 0 0;0 1 0;-1 0 0;0 0 0];
        
    else
        % make it a column vec. if necess.
        if n==length(IN)
            IN=IN';
        end
        
        dm3dIN = [eye(3);IN'/theta];
        omega = IN/theta;
        dm2dm3 = [eye(3)/theta -IN/theta^2; zeros(1,3) 1];
        
        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1-cos(theta);
        omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        A = omega*omega';
        
        dm1dm2 = zeros(21,4);
        dm1dm2(1,4) = -sin(theta);
        dm1dm2(2,4) = cos(theta);
        dm1dm2(3,4) = sin(theta);
        dm1dm2(4:12,1:3) = [0 0 0 0 0 1 0 -1 0;0 0 -1 0 0 0 1 0 0;0 1 0 -1 0 0 0 0 0]';
        
        w1 = omega(1);
        w2 = omega(2);
        w3 = omega(3);
        
        dm1dm2(13:21,1) = [2*w1;w2;w3;w2;0;0;w3;0;0];
        dm1dm2(13: 21,2) = [0;w1;0;w1;2*w2;w3;0;w3;0];
        dm1dm2(13:21,3) = [0;0;w1;0;0;w2;w1;w2;2*w3];
        
        R = eye(3)*alpha + omegav*beta + A*gamma;
        
        dRdm1 = zeros(9,21);
        dRdm1([1 5 9],1) = ones(3,1);
        dRdm1(:,2) = omegav(:);
        dRdm1(:,4:12) = beta*eye(9);
        dRdm1(:,3) = A(:);
        dRdm1(:,13:21) = gamma*eye(9);
        dRdIN = dRdm1 * dm1dm2 * dm2dm3 * dm3dIN;
        
        
    end
    OUT = R;
    dOUT = dRdIN;
    
% input is rotational matrix
elseif ((m==n) && (m==3) && (norm(IN' * IN - eye(3)) < bigeps) && (abs(det(IN)-1) < bigeps))
    
    R = IN;
    
    % project the rotation matrix to SO(3)
    [U,~,V] = svd(R);
    R = U*V';
    tr = (trace(R)-1)/2;
    dtrdR = [1 0 0 0 1 0 0 0 1]/2;
    theta = real(acos(tr));
    
    if sin(theta) >= 1e-4
        
        dthetadtr = -1/sqrt(1-tr^2);
        dthetadR = dthetadtr * dtrdR;
        vth = 1/(2*sin(theta));
        dvthdtheta = -vth*cos(theta)/sin(theta);
        dvar1dtheta = [dvthdtheta;1];
        dvar1dR =  dvar1dtheta * dthetadR;
        om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
        dom1dR = [0 0 0 0 0 1 0 -1 0;0 0 -1 0 0 0 1 0 0;0 1 0 -1 0 0 0 0 0];
        dvardR = [dom1dR;dvar1dR];
        om = vth*om1;
        domdvar = [vth*eye(3) om1 zeros(3,1)];
        dthetadvar = [0 0 0 0 1];
        dvar2dvar = [domdvar;dthetadvar];
        OUT = om*theta;
        domegadvar2 = [theta*eye(3) om];
        dOUT = domegadvar2 * dvar2dvar * dvardR;
        
    else
        % case norm(v_rot)=0;
        if tr > 0;
            
            OUT = [0 0 0]';
            dOUT = [0 0 0 0 0 1/2 0 -1/2 0; 0 0 -1/2 0 0 0 1/2 0 0;0 1/2 0 -1/2 0 0 0 0 0];
            
        % case norm(v_rot)=pi;   
        else
            
            % define hashtable
            hash_tab = [0; -1; -3; -9; 9; 3; 1; 13; 5; -7; -11];
            Smat = [1,1,1; 1,0,-1; 0,1,-1; 1,-1,0; 1,1,0; 0,1,1; 1,0,1; 1,1,1; 1,1,-1;1,-1,-1; 1,-1,1];
            
            M = (R+eye(3,3))/2;
            uabs = sqrt(M(1,1));
            vabs = sqrt(M(2,2));
            wabs = sqrt(M(3,3));
            
            mvec = [M(1,2), M(2,3), M(1,3)];
            syn  = ((mvec > 1e-4) - (mvec < -1e-4));
            hash = syn * [9; 3; 1];
            idx = find(hash == hash_tab);
            svec = Smat(idx,:)';
            
            OUT = theta * [uabs; vabs; wabs] .* svec;
            
        end
    end
end

