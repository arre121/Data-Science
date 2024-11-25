function [m] = RK3(diff,start_val,T,N)
t(1)=0;
h=T/N;
m(:,1)=start_val;
    for i=1:N
        %Update time
        t(i+1)=t(i)+h;
        %Update m
        k1=diff(t(i),m(:,i));
        k2=diff(t(i)+h,m(:,i)+h*k1);
        k3=diff(t(i)+h/2,m(:,i)+h*k1/4+h*k2/4);
        m(:,i+1)=m(:,i)+h/6*(k1+k2+4*k3);
    end

end

