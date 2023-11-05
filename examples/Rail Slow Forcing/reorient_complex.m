function VC = reorient_complex(M,M0)
V_int = [M(:,1),M(:,2),M(:,3)]
Theta = [];
V_comp = [M0(:,1),M0(:,2),M0(:,3)];
VC = [];
for indf = 1:max(size(V_int(1,:)))
Df = V_int(:,indf)./V_comp(:,indf);
TF = isinf(Df);
Df(find(TF)) = 0;
TF = isnan(Df);
Df(find(TF)) = 0;
TF = find(Df);
d = Df(TF(1));
[theta,rho] = cart2pol(real(d),imag(d));
%if abs(theta) >= pi/2
Vchange = exp(-1i*theta).*V_int(:,indf);
%else
%Vchange = V_int(:,indf);
%end    
VC =[VC,Vchange];
end    
VC =[VC,conj(Vchange)];
end