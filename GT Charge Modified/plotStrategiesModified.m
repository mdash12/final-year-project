clear;
axis([0 500 0 500]);
hold on;
plot(250,250,'bs','LineWidth',3);



noOfNodes=textread('Put_NO_OF_NODES.txt','%d');
disp(noOfNodes);


node_X=textread('Put_NODE_X.txt','%f',noOfNodes);
node_Y=textread('Put_NODE_Y.txt','%f',noOfNodes );

plot(node_X,node_Y,'r*','LineWidth',3); 


no_of_chargers=textread('Put_NO_OF_CHARGERS.txt','%d');
disp(no_of_chargers);
charger_X=textread('Put_CHARGER_X.txt','%f',no_of_chargers);
charger_Y=textread('Put_CHARGER_Y.txt','%f',no_of_chargers);

plot(charger_X,charger_Y,'gV','LineWidth',3); 


no_of_games=textread('Put_NO_OF_GAMES.txt','%d');
disp(no_of_games);
x= textread('Put_STRATEGY.txt');
%disp(x);
no_of_games=no_of_games+1;
c=1/no_of_chargers;
%disp(c);

for i=1:4
   temp=x(((i-1)*no_of_games+1):(i*no_of_games),:);
   k=(i-1)*c;
   %disp(k);
   %disp(k);
   for j=1:no_of_games-1
       
        p1 = [temp(j,1) temp(j,2)];       % First Point
        p2 = [temp(j+1,1) temp(j+1,2)];       % Second Point
        dp = p2-p1;                         % Difference
        quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',1.5,'Color',[k*k 1-k k])
        
       
       aa=int2str(uint32(temp(j+1,3)));
       xx=(temp(j,1)+temp(j+1,1))/2;
       yy=(temp(j,2)+temp(j+1,2))/2;
       dis =sqrt((temp(j,1)-temp(j+1,1))*(temp(j,1)-temp(j+1,1)) + ((temp(j,2)-temp(j+1,2))*(temp(j,2)-temp(j+1,2))));
       strDis=num2str(dis);
       %disp(dis);
       %g=int2str(uint32(j))
       text(xx,yy,int2str(uint32(j)));
       %text(xx,yy,aa);
       hold on;
       
   end
end
   