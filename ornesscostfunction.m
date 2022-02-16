function hh=costfunction(x)
%  global grst746
%  global cmst746 msst746 prst746 trst746rt trst746V7

x=x';

for i=1:4
         orness(i)=-x(i)*log(x(i));
end
z=sum(orness(:));

hh=z;

