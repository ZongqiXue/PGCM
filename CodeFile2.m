clc; clear all;
load Dataset2.mat
theta = [5 5]; lob = [1e-1 1e-1]; upb = [20 20];
X = gridsamp([140.994 29.591;141.117 29.743], 80);
X1 = reshape(X(:,1),80,80); X2 = reshape(X(:,2),80,80);
YX=zeros(length(X),35);
for i=1:35
    data=[bsgi(1:311,2:3),bsfit(1:311,i)];
    eval(['[dmodel',num2str(i),' perf',num2str(i),'] = dacefit(data(:,1:2), data(:,3), @regpoly0, @corrgauss, theta, lob, upb);']);
    eval(['[YX(:,i) MSE] = predictor(X, dmodel',num2str(i),');']);
end

for t=1:12
    plotkriging(bsfit,bsgi,YX,X1,X2,t);

end

