clear all;clc;
format compact;
fprintf("\n\n===============FOR FUTURE CRUDE OIL===============\n\n")
%loading data for Future Crude Oil
[priceCO, dateCO, wholeCO]= xlsread('Assignment02 Data.xlsx',"Future Crude Oil");
%return of Future Crude Oil
rawreturnCO = price2ret(priceCO);
%split into two parts
returnCO = rawreturnCO(1:length(rawreturnCO)-20); %for modeling
returnCO2 = rawreturnCO(length(rawreturnCO)-19 :length(rawreturnCO)); %for estimating forecast

%trials and error: find the suitable model for Future Crude Oil
highest_aic=0;
highest_aic_name="";
highest_bic=0;
highest_bic_name="";
%A = ["ARIMA(p,q)" "AIC" "BIC"];
A = [];
B = [];
C = [];
ctr = 0;
for p = 0:5
for q = 0:5
% switch between normal dist & t-dist
for N_T_switch = 1:2
model = arima(p,0,q);
if N_T_switch == 2
model.Distribution = "t" ;
fprintf("\n\n ============== The below is t-distribution ============\n\n");
end
fprintf("arima(%i,0,%i)", p,q);
ESTmodelN = estimate(model,returnCO);
[res0N,v0N,logL0N] = infer(ESTmodelN,returnCO);
[aic, bic] = aicbic(logL0N, p+q+1);
% t-dist has one more parameter to calculate
if N_T_switch==2
[aic, bic] = aicbic(logL0N,p+q+2);
end
stdres0N = res0N./sqrt(v0N);
[h,prob]=lbqtest(stdres0N,'Lags',[5,10,15])
disp(aic)
disp(bic)
% if there is no significant correlation, compare aicbic
if h(1) == 0 %eliminate the model(s) that cannot completely
pass Ljung-Box Test
name = "arima(" + int2str(p) + ",0," + int2str(q) + ")";
if N_T_switch ==2
name = name + "t-dist";
end
A = [A ; name];
B = [B ; aic];
C = [C ; bic];
end
end
end
end
[B,I1] = sort(B);
A = A(I1);
fprintf("AIC performance of Models for Future Crude Oil: Smallest to
Highest \n")
[A sort(B)]
[C,I2] = sort(C);
A = A(I2);
fprintf("BIC performance of Models for Future Crude Oil: Smallest to
Highest \n")
[A sort(C)]
%trials and error: find the suitable model for Future Crude Oil
fprintf("\n\n===============FOR SNP500===============\n\n")
%loading data for S&P500
[priceSP, dateSP, wholeSP]= xlsread('Assignment02 Data.xlsx',"S&P500");
%return of S&P500
rawreturnSP = price2ret(priceSP);
%split into two parts
returnSP = rawreturnSP(1:length(rawreturnSP)-
20); %for modeling
returnSP2 = rawreturnSP(length(rawreturnSP)-19 :
length(rawreturnSP)); %for estimating forecast
highest_aic=0;
highest_aic_name="";
highest_bic=0;
highest_bic_name="";
A = []; %A = ["ARIMA(p,q)"]
B = []; %B = [aic_value]
C = []; %B = [bic_value]
ctr = 0;
for p = 0:5
  for q = 0:5
    % switch between normal dist & t-dist
    for N_T_switch = 1:2
      model = arima(p,0,q);
      if N_T_switch == 2
      model.Distribution = "t" ;
      fprintf("\n\n ============== The below is t-distribution
      ============\n\n");
    end
    fprintf("arima(%i,0,%i)", p,q);
    ESTmodelN = estimate(model,returnSP);
    [res0N,v0N,logL0N] = infer(ESTmodelN,returnSP);
    [aic, bic] = aicbic(logL0N, p+q+1);
    % t-dist has one more parameter to calculate
    if N_T_switch==2
      [aic, bic] = aicbic(logL0N,p+q+2);
      end
      stdres0N = res0N./sqrt(v0N);
      [h,prob]=lbqtest(stdres0N,'Lags',[5,10,15])
      disp(aic)
      disp(bic)
      % if there is no significant correlation, compare aicbic
      if h(1) == 0 %clean off the model(s) that cannot pass Ljung-Box
        Test
        name = "arima(" + int2str(p) + ",0," + int2str(q) + ")";
        if N_T_switch ==2
          name = name + "t-dist";
        end
        A = [A ; name];
        B = [B ; aic];
        C = [C ; bic];
      end
    end
  end
end
[B,I1] = sort(B);
A = A(I1);
fprintf("AIC performance of Models for SNP500: Smallest to Highest \n")
[A sort(B)]
[C,I2] = sort(C);
A = A(I2);
fprintf("BIC performance of Models for SNP500: Smallest to Highest \n")
[A sort(C)]
%forecasting
%for Future Crude Oil ARMA(4,3)
modelCO = arima(4,0,3)
ESTmodelCO = estimate(modelCO,returnCO)
[E0CO,V0CO] = infer(ESTmodelCO, returnCO);
%95% interval forecasts
[YCO,YMSECO,VCO] = forecast(ESTmodelCO,20,returnCO,"E0",E0CO,"V0",V0CO);
upperCO = YCO + 1.96*sqrt(YMSECO);
lowerCO = YCO - 1.96*sqrt(YMSECO);
figure
plot(returnCO)
hold on
plot(length(returnCO)+1:length(returnCO)+20,YCO,"r","LineWidth",2)
plot(length(returnCO)+1:length(returnCO)+20,[upperCO,lowerCO],"k--")
xlim([0,length(returnCO)+20])
title("Forecasted Returns of Future Crude Oil")
hold off
%mse and mae
mseCO = mse(returnCO2-YCO)
maeCO = mae(returnCO2-YCO)
%for S&P500 ARMA(1,1)
modelSP = arima(1,0,1)
ESTmodelSP = estimate(modelSP,returnSP)
[E0SP,V0SP] = infer(ESTmodelSP, returnSP);
%95% interval forecasts
[YSP,YMSESP,VSP] = forecast(ESTmodelSP,20,returnSP,"E0",E0SP,"V0",V0SP);
upperSP = YSP + 1.96*sqrt(YMSESP);
lowerSP = YSP - 1.96*sqrt(YMSESP);
figure
plot(returnSP)
hold on
plot(length(returnSP)+1:length(returnSP)+20,YSP,"r","LineWidth",2)
plot(length(returnSP)+1:length(returnSP)+20,[upperSP,lowerSP],"k--")
xlim([0,length(returnSP)+20])
title("Forecasted Returns of S&P500")
hold off
%mse and mae
mseSP = mse(returnSP2-YSP)
maeSP = mae(returnSP2-YSP)
