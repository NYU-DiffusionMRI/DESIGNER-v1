function dirs = radialsampling(dir, n)

% compute Equator Points

dt = 2*pi/n;
theta = 0:dt:(2*pi-dt);

dirs = [cos(theta)', sin(theta)', 0*theta']';

v = [-dir(2), dir(1), 0];
s = sqrt(sum(v.^2));
c = dir(3);
V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + V + V*V * (1-c)/s^2;

dirs = R*dirs;

end