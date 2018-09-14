clear
clc

% Gone for 24 days ... (3.42 weeks)
wks_vacation = 24/7;
hrs_vacation = 40 * wks_vacation;

%%% 40 hours/wk for May, June, July (3 months) minus vacation time
weeksSummer = 3*30.6/7;
hoursWorked_40 = 40 * weeksSummer - hrs_vacation


%%% 30 hours/wk 
hoursWorked_30 = 30 * weeksSummer