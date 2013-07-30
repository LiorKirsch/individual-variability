function [area,n_true_pos,n_fals_pos,sens,spec,acc,threshes]=G_roc(true,predicted)
%
% ROC 
% Purpose:calculate roc curve
%
% Usage: [AUC,TP,FP,SENS,SPEC,ACC,THRESH]=ROC(TRUE,PREDICTED)
%
% Inputs: 
%    predicted  - real values in [0,1]
%    true       - true labels (binary) {0,1}  
%
%    All values of 'true' smaller than 0 are set to 0 
%    (this effectively transforms {-1,1}==>{0,1} ) 
%       
% Outputs: 
%    AUC      - Area under teh curve (AUC)
%    TP       - No of true positives for each threshold value
%    FP       - No of true positives for each threshold value  
%    SENS     - Sensitivity: T_POS / (T_POS + F_NEG)
%    SPEC     - Specificity: T_NEG / (T_NEG + F_POS)
%    ACC      - Accuracy for each threshold value
%    THRESHES - The threshold that achieves the best accuracy
%
% To plot the ROC curve use: 
%    plot(tp,1-fp)
%    xlabel('True positives');
%    ylabel('1 - False positives');
% 
  
  n = length(true);
  if(length(predicted)~=n)
    error('Illegal input: Inconsistent input to eoc.m: %d  vs %d\n',n,length(predicted));
    n_false_pos=0;
    n_true_pos =0;
    area = 0;
    return;
  end
  % Make sure this is a column vector
  predicted = predicted(:);
  true      = true(:);  
  
  true(find(true<0))=0;
  
  n_true1 = length(find(true));
  n_true0 = length(find(~true));
  if(n_true0==0 || n_true1 ==0)
    error('G_roc:true_same','Illegal input: True-labels are all the same');    
  end  
  
  
  threshes = [sort(predicted-eps*2);max(predicted)+eps*2];
  n_thr =length(threshes);
  
  binary_pred = zeros(1,n);
  n_true_pos = zeros(1,n_thr);
  n_fals_pos = zeros(1,n_thr);
  n_true_neg = zeros(1,n_thr);
  n_fals_neg = zeros(1,n_thr);
  sens = zeros(1,n_thr);
  spec = zeros(1,n_thr);
  
  for i = 1:n_thr
    thresh  = threshes(i);
    binary_pred = (predicted>thresh);
    n_true_pos(i) = length(find(binary_pred==1 & true==1));
    n_fals_pos(i) = length(find(binary_pred==1 & true==0));
    n_true_neg(i) = length(find(binary_pred==0 & true==0));
    n_fals_neg(i) = length(find(binary_pred==0 & true==1));
    sens(i) = n_true_pos(i) / (n_true_pos(i)+ ...
					     n_fals_neg(i));
    spec(i) = n_true_neg(i) / (n_true_neg(i)+ ...
					     n_fals_pos(i));   
  end
  
  acc = (n_true_pos + n_true_neg) / n;
  

  area = 0;
  for(i = 2:n_thr);
    % area=area+((1-spec(i))+(1-spec(i-1)))/2*(sens(i)-sens(i-1));
    x = (1-spec(i))-(1-spec(i-1));
    area=area+((1-spec(i))-(1-spec(i-1)))*(sens(i)+sens(i-1))/2;    	 
  end

  area = abs(area);

  if(isnan(area)), keyboard; error; end;
  
  if(nargout<1)
    plot(1-spec,sens);
    xlabel('1 - specifcity');
    ylabel('sensitivity');    
  end

  
  
  return
