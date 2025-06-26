function [dfctr,delta] = controlDispFctr(dfctr,delta)

dfctr  = dfctr - delta;
delta =  delta * 0.1;
dfctr =  dfctr + delta;


end
