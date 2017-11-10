%HW6

% Problem 1. Curve fitting. 
% Part 1. Take 10 x values on the interval 0 < x <= 10 and then create y
% values from the x values by plugging the x values into a third order
% polynomial of your choice. Add random noise to the data by choosing a random number
% in the interval [-D, D]. Start with D = 1. Plot your data in the x-y plane.

X = sort(10*rand(1,10));

Y= X.^3+2*X+1;
noise=randn(1,10);
Y1=Y+noise;
figure
plot(X,Y1)
% Part 2. Fit your data with polynomials from order 1 to 9. Plot the fitted
% polynomials together with the data. 
figure
plot(X,Y1)
hold on
for i = 1:9
    [coeff, s] = polyfit(X,Y1,i);
    yfit = polyval(coeff,X);
    plot(X,yfit);
    hold on
end
% Part 3. On a separate plot, plot the R^2 and adjusted R^2 as a function
% of the order of the polynomial. 
R=[];
R1=[];
for i=1:9
    [coeff, s]=polyfit(X,Y1,i);
    yfit=polyval(coeff,X);
    Yr=Y1-yfit;
    ssresid=sum(Yr.^2);
    sstotal=(length(Y1)-1)*1;
    R(i)=1-ssresid/sstotal;
    R1(i)=1-ssresid/sstotal*(length(Y1)-1)/(length(Y1)-length(coeff));
end
figure
plot(1:9,R);
hold on
plot(1:9,R1);


% Part 4. Repeat parts 1 - 3 for D = 10 and D = 1000. Comment on the
% results. 

%D=10
X = sort(10*rand(1,10));

Y= X.^3+2*X+1;
noise=10*randn(1,10);
Y1=Y+noise;
figure
plot(X,Y1)
hold on;
R=[];
R1=[];
for i=1:9
    [coeff, s]=polyfit(X,Y1,i);
    yfit=polyval(coeff,X);
    Yr=Y1-yfit;
    s=sum(Yr.^2);
    s1=(length(Y1)-1)*10;
    R(i)=1-s/s1;
    R1(i)=1-s/s1*(length(Y1)-1)/(length(Y1)-length(coeff));
    hold on
    plot(X,yfit);
end
figure(2)
plot(1:9,R);
hold on;
plot(1:9,R1);

%D=100
X = sort(10*rand(1,10));

Y= X.^3+2*X+1;
noise=100*randn(1,10);
Y1=Y+noise;
figure
plot(X,Y1)
hold on;
R=[];
R1=[];
for i=1:9
    [coeff, s]=polyfit(X,Y1,i);
    yfit=polyval(coeff,X);
    Yr=Y1-yfit;
    s=sum(Yr.^2);
    s1=(length(Y1)-1)*100;
    R(i)=1-s/s1;
    R1(i)=1-s/s1*(length(Y1)-1)/(length(Y1)-length(coeff));
    hold on
    plot(X,yfit);
end

figure(2);
plot(1:9,R);
hold on;
plot(1:9,R1);

% The increasing noise in the model makes the lower orders of polynomial
% fits worse than the higher orders (3 and above). The larger the noise is
% ,the worse the fitting curve of the lower orders is. Generally, higher
% orders polynomial fits the data better.

% Part 5. Now repeat parts 1-3 but take 100 x values on the interval 0 < x <=
% 10. Comment on the results. 
%D=10
X = sort(10*rand(1,100));

Y= X.^3+2*X+1;
noise=randn(1,100);
Y1=Y+noise;
figure
plot(X,Y1)
hold on;
R=[];
R1=[];
for i=1:9
    [coeff, s]=polyfit(X,Y1,i);
    yfit=polyval(coeff,X);
    Yr=Y1-yfit;
    s=sum(Yr.^2);
    s1=(length(Y1)-1);
    R(i)=1-s/s1;
    R1(i)=1-s/s1*(length(Y1)-1)/(length(Y1)-length(coeff));
    hold on
    plot(X,yfit);
end

figure(2)
plot(1:9,R);
hold on
plot(1:9,R1);

%100x values curve is smoother than the 10x values. The adjusted Rsquare is also closer to the unadjusted.

% Problem 2. Basic statistics. 
% Part 1. Consider two different distributions - Gaussian numbers with a mean of
% 0 and variance 1 and Gaussian numbers with a mean of 1 and variance 1.
% (1) Make a plot of the average p-value for the t-test comparing N random
% numbers chosen from each of these two distributions as a function of N.

function plotgaussian(N)
d=[];
for i=1:N
        a=randn(N);
        b=randn(N)+1;
        [is_sig,pval]=ttest2(a,b);
    d(i)=mean(pval);
end
xlab=1:N;
figure
plot(xlab,d)
end

plotgaussian(5)
% Part 2. Now keep the first distribution the same, but vary the mean of
% the second distribution between 0 and 10 with the same variance and
% repeat part one. Make a plot of all of these different curves on the same
% set of axes. What is special about the case where the mean of the second
% distribution is 0? 
N=100;
d=[]
for M=0:10
    for i=1:N
            a=normrnd(0,1,[1,N]);
            b=normrnd(M,1,[1,N]);
            [~,d(M+1,i)]=ttest2(a,b);
    end
end

for M=1:10
xlab=1:N;
figure
    plot(xlab,d(M+1,:));
    hold on;
end

% when the mean of the second distriution is 0, which the same as the first
% distribution, the mean of pval is significantly higher than others.So we
% can only see that curve given how the axis are set. This means we cannot
% reject the null when the second distribution's mean is 0. It also means
% the 2 distributions are the same.

% Part 3. Now keep the means of the two distributions at 0 and 1 as in part
% 1, but vary the variance of both distributions simultaneiously between 0.1 and 10 and plot the 
% p-values vs the number of numbers drawn as before. Comment on your results.  

N = 100;
p_val_m2 = [];
for var = 0.1:1:10.1
    for ii = 1:N
        A = normrnd(0,var,[1,N]);
        B = normrnd(1,var,[1,N]);
        [~,p_val_m2((var-0.1)+1,ii)] = ttest2(A,B);
    end
end
figure
title('p-values vs. N');

for i = 1:11
    plot(1:N,p_val_m2(i,:))
    hold on
end

% when the variance is small, the p values are small as well. This is
% because the distribution of the second distribution are all around the
% mean. Hence it has a greater chance to reject the null when the variance
% is very small. As the variance increases, the p-value increases. because
% the higher the variance is, the more chance the number in the second
% distribution will overlap the first one, and the higher possibility that
% we cannot reject the null.
