function Vp = U2plateau(V, Sheffield_sequence)
% V could be based on Sheffield voltage measurement sequence
% or other sequence like measuring followed by injecting electrodes
%reshape sinusoidal wave form (Vs) to plateau shaped wave form (Vp)

% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3

for j=1:size(V,2)
    if Sheffield_sequence
        Vs = V(:,j);
        k=0;
        for m=1:16
            if m <= 2 || m ==16
                Vm((m-1)*13+1:m*13)= Vs((m-1)*13+1:m*13);
            else
                k =1+k;
                for count=1:k
                    Vm(m*13+count-k)=Vs((m-1)*13+count);
                end
                Vm((m-1)*13+1:m*13-k)= Vs((m-1)*13+1+k:m*13);
            end
        end
    else
        Vm = V(:,j);
    end
    
    for m = 1:13
        index1((m-1)*16+1:m*16) = m:13:13*16+m-1;
    end
    Vp(:,j) = Vm(index1);
end
