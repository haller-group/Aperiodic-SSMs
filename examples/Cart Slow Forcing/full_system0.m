function dy = full_system0(t,y,epsilon,F)
m1 = 1;
m2 = 1;
Mf = 2;
kf = 1;
k = 1;
c = 0.3;
cf = 0.3;
gamma = 0.5;
MTT = m1 + m2 + Mf;
K = [-2*k - kf*((m1 + m2)/MTT)^2, -k - kf*m2*(m1 + m2)/MTT^2, kf*(m1+m2)/MTT; -k - kf*m2*(m1 + m2)/MTT^2, -2*k - kf*m2^2/(MTT^2), kf*m2/MTT; (m1 + m2)/MTT *kf, m2/MTT*kf, -kf];
C = [-c - cf*((m1 + m2)/MTT)^2, -c - cf*m2*(m1 + m2)/MTT^2, cf*(m1+m2)/MTT; -c - cf*m2*(m1 + m2)/MTT^2, -2*c - cf*m2^2/(MTT^2), cf*m2/MTT; (m1 + m2)/MTT *cf, m2/MTT*cf, -cf];
M =  [Mf*(m1 + m2)/MTT, m2*Mf/MTT, 0; m2*Mf/MTT, m2*(m1 + Mf)/MTT, 0; 0,0, MTT];


A = [zeros(3),eye(3);inv(M)*K,inv(M)*C];
Linpart = A*y;
dy=Linpart + [0;0;0;-gamma*(m1+Mf)/(m1*Mf)*y(1)^3;gamma/m1*y(1)^3;0];
end