%Assume that we are given a position vector, in which the collection times
 %are (approximately) equidistantly spaced, with distance h. 
 %The higher the order of the approximation, the more data we will lose
 %from the ends. However, this is generally negligible. the 2nd order
 %derivative loses 2 total data points, 1 from each end.
 function y=velocity2(x,P)
 h=mean(x(2:end)-x(1:end-1));

 y=(P(3:end)-P(1:end-2))./(2*h);
 