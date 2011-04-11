function RERi = recip_err(D1, D2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:recip_err.m
% RER --> Reciprocity Error calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (C) 2011 Mamatjan Yasheng. License: GPL v2 or v3

Max1=find(D1==max(D1));
Max2=find(D2==max(D2));
MD = abs(Max1-Max2);
L=length(D1);
D22(1:L-MD)=D2(MD+1:L);
D22(L-MD+1:L)=D2(1:MD);

figure
subplot(2,1,1)
% plot to see if the figures match or not
plot(D1,'DisplayName','D1','YDataSource','D1');
hold all;
plot(D22,'DisplayName','D22','YDataSource','D22');
hold all;hold off;figure(gcf);
axis tight

tmp_sum = 0;
for i =1:L
    tmp(i) = abs(D1(i)-D22(i))/D1(i);
    tmp_sum = tmp_sum + (tmp(i))^2;
end
RER_sum = sqrt(tmp_sum/L)

IDX = [18:2:416,2:2:16];
RERi= (1-tmp(IDX))*100;

subplot(2,1,2)
plot(RERi); axis tight
xlabel('Channel number');
ylabel('Reciprocity accuracy');
end
