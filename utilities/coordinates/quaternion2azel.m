function [az,el,rot] = quaternion2azel(qx,qy,qz,qw)
     
    n=numel(qx);
    az = zeros(1,n);
    el = zeros(1,n);
    rot = zeros(1,n);
    
    for k=1:n
        q = quaternion(qx(k), qy(k), qz(k), qw(k));
        qRot = quaternion(0, 0, 0, 1);
        q=q*qRot;
        a = q.EulerAngles('yzx');
        az(k) = a(1) * -180.0 / pi;  % AZ
        el(k) = a(2) * -180.0 / pi;  % EL
        rot(k) = a(3) * 180.0 / pi;  % ROT
    end
end

