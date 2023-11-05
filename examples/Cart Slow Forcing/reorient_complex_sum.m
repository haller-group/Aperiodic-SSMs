function VC = reorient_complex_sum(M)
V_int = M(:,1:2:end);
VC = [];
for indf = 1:max(size(V_int(1,:)))
Df = sum(V_int(:,indf));
if (real(Df) > 0 && imag(Df)<0) 
    Vchange = 1i*V_int(:,indf);
elseif (real(Df) < 0 && imag(Df)<0) 
    Vchange = -V_int(:,indf);
elseif (real(Df) < 0 && imag(Df)>0) 
    Vchange = -1i*V_int(:,indf);
else
    Vchange = V_int(:,indf);
end
VC =[VC,Vchange,conj(Vchange)];
end    

end