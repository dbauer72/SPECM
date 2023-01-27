A = [1,0;0,.5];
K = [1,0;0,1];C= [1,0;0,1];
Abar = A-K*C;

beta = [0,1]';
alpha = [0;1/3];

T= 50
M = 1000
[alphahat,alphahat_1,alphahat_2,alphahat_3,alphahat_v,alpha] = simu_Calpha(A,K,C,beta,T,M);

ahat_50 = [alphahat(:,1),alphahat_1(:,1),alphahat_2(:,1),alphahat_3(:,1),alphahat_v(:,1)];
that_50 = [alphahat(:,3),alphahat_1(:,3),alphahat_2(:,3),alphahat_3(:,3),alphahat_v(:,3)];

T= 100
M = 1000
[alphahat,alphahat_1,alphahat_2,alphahat_3,alphahat_v,alpha] = simu_Calpha(A,K,C,beta,T,M);

ahat_100 = [alphahat(:,1),alphahat_1(:,1),alphahat_2(:,1),alphahat_3(:,1),alphahat_v(:,1)];
that_100 = [alphahat(:,3),alphahat_1(:,3),alphahat_2(:,3),alphahat_3(:,3),alphahat_v(:,3)];


T= 200
M = 1000
[alphahat,alphahat_1,alphahat_2,alphahat_3,alphahat_v,alpha] = simu_Calpha(A,K,C,beta,T,M);

ahat_200 = [alphahat(:,1),alphahat_1(:,1),alphahat_2(:,1),alphahat_3(:,1),alphahat_v(:,1)];
that_200 = [alphahat(:,3),alphahat_1(:,3),alphahat_2(:,3),alphahat_3(:,3),alphahat_v(:,3)];

T= 500
M = 1000
[alphahat,alphahat_1,alphahat_2,alphahat_3,alphahat_v,alpha] = simu_Calpha(A,K,C,beta,T,M);

ahat_500 = [alphahat(:,1),alphahat_1(:,1),alphahat_2(:,1),alphahat_3(:,1),alphahat_v(:,1)];
that_500 = [alphahat(:,3),alphahat_1(:,3),alphahat_2(:,3),alphahat_3(:,3),alphahat_v(:,3)];

save ahat_simu

[nanmean(ahat_50);nanmean(ahat_100);nanmean(ahat_200);nanmean(ahat_500)]
[mean(abs(ahat_50).^2);mean(abs(ahat_100).^2);nanmean(abs(ahat_200).^2);nanmean(abs(ahat_500).^2)]
[nanmean(that_50);nanmean(that_100);nanmean(that_200);nanmean(that_500)]
