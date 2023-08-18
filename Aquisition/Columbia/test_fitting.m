% Construct the dark noise and hardware stability regressors

scan_ref=load('RefSpectrum-Patient1_20221221T101545_scan7.mat');
scan_dark=load('RefSpectrum-Patient1_20221221T111940_scan8.mat');
scan_dark2=load('RefSpectrum-Patient1_20221221T122137_scan9.mat');

scan_dark.raw_data=nirs.math.tddr(scan_dark.raw_data',5)';
scan_dark2.raw_data=nirs.math.tddr(scan_dark2.raw_data',5)';
scan_ref.raw_data=nirs.math.tddr(scan_ref.raw_data',5)';

scan_dark.raw_data=nirs.math.tddr(scan_dark.raw_data,5);
scan_dark2.raw_data=nirs.math.tddr(scan_dark2.raw_data,5);
scan_ref.raw_data=nirs.math.tddr(scan_ref.raw_data,5);

ref.dark(:,1) = median(scan_dark.raw_data,1);
ref.dark(:,2) = median(scan_dark2.raw_data,1);
ref.dark(:,3)=1;
ref.light(:,1)= median(scan_ref.raw_data,1);

fluct1=scan_dark.raw_data-ones(size(scan_dark.raw_data,1),1)*ref.dark(:,1)';
fluct2=scan_dark2.raw_data-ones(size(scan_dark2.raw_data,1),1)*ref.dark(:,2)';
fluct3=scan_ref.raw_data-ones(size(scan_ref.raw_data,1),1)*ref.light(:,1)';

[~,s1,v1]=nirs.math.mysvd(cov(fluct1));
[~,s2,v2]=nirs.math.mysvd(cov(fluct2));
[~,s3,v3]=nirs.math.mysvd(cov(fluct3));

nSV1=min(find(cumsum(diag(s1))/sum(diag(s1))>.8));
nSV2=min(find(cumsum(diag(s2))/sum(diag(s2))>.8));
nSV3=min(find(cumsum(diag(s3))/sum(diag(s3))>.8));


ref.PCA = orth([ref.dark v1(:,1:nSV1) v2(:,1:nSV2) v3(:,1:nSV3)]);

refd=ted_BH.raw_data;
refd=nirs.math.tddr(refd',5)';
refd=nirs.math.tddr(refd,5);
%I = A*Io*exp(-mua)
% mua = log(A)+log(I0)-log(I)
refd = log(ones(size(refd,1),1)*abs(ref.light)')-log(refd)+log(3000);

mrefd=median(refd,1);
refd=(refd'-PCA*inv(PCA'*PCA)*PCA'*refd')';
refd=refd+ones(size(refd,1),1)*mrefd;


