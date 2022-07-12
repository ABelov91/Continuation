figure;
hold on;
plot(t,u,'-k','MarkerFaceColor','k','Linewidth',1)
plot(t,v,'-r','MarkerFaceColor','r','Linewidth',1)
plot(t,w,'-ob','MarkerFaceColor','b','Linewidth',1)
plot(t,rh,'-g','MarkerFaceColor','b','Linewidth',1)
plot(t,order_u,'-m')
plot(t,flag,':k')
axis([0 TIME -5 15]);
legend('u', 'v', 'w', 'f', 'order u', 'flag')