function r=eno3c(um2,um1,u,up1,up2,~)
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


end
