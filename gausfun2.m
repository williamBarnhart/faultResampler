function res=gausfun2(in)
global x2 b y2
id2=find(x2==0);

if isempty(id2)
    id2=find(x2==min(x2));
end
%a=in(1);
c=in(1);
%res=norm(a*exp(-((x2-b)/c).^2)-y2);

res=y2(id2)*exp(-((x2-b)/c).^2)-y2;

