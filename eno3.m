% To generate regression data using eno3
function ul=eno3(um2,um1,u,up1,up2,~)
c= [1/3,-7/6,11/6;-1/6,5/6,1/3;1/3,5/6,-1/6];
state=[um2, um1, u, up1, up2];
rec_order=3;
ul=0;
dd=zeros(rec_order,2*rec_order-1);
for j=1:2*rec_order-1
    dd(1,j) = state(j);
end
for i=2:rec_order
    for j=1:2*rec_order- i
        dd(i,j) = dd(i-1,j+1) - dd(i-1,j);
    end
end
r=0;
s = rec_order-1;
for i=2:rec_order
    if(abs(dd(i,s-r+1)) > abs(dd(i,s-r)))
        r=r+1;
    end
end
for i=1:rec_order
    ul =ul+ c(s-r+1,i)*state(s-r+i);
end

end
