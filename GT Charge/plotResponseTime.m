clear;

no_of_games=textread('Put_NO_OF_GAMES.txt','%d');
% disp(no_of_games);

no_of_nodes=textread('Put_NO_OF_NODES.txt','%d');
% disp(no_of_nodes);

axis([0 no_of_nodes 0 no_of_games]);
hold on;


reqGame=textread('Put_REQ_GAME.txt');
servedGame=textread('Put_SERVED_GAME.txt');


reqGameModified=textread('F:\\project8\\GT Charge Modified\\Put_REQ_GAME.txt');
servedGameModified=textread('F:\\project8\\GT Charge Modified\\Put_SERVED_GAME.txt');

%GTCharge
for i=1:no_of_nodes
   temp1=reqGame(i:i,:);
   temp2=servedGame(i:i,:);
   
   size1=temp1(1);
   size2=temp2(1);
   
   k=(i-1)*1/no_of_nodes;
   
   line([i i],[0 no_of_games],'Color','y');
   hold on;
   
   for j=2:size1+1
       plot(i,temp1(j),'r*','LineWidth',0.5);
       hold on; 
   end
   
   hold on;
   
   for j=2:size2+1
       if(temp2(j)==temp1(j))
            plot(i,temp2(j),'b*','LineWidth',0.5);
       else
            plot(i,temp2(j),'g*','LineWidth',0.5);
       end
       hold on; 
   end    
end


%GTCharge Modified
% for i=1:no_of_nodes
%    temp1=reqGameModified(i:i,:);
%    temp2=servedGameModified(i:i,:);
%    
%    size1=temp1(1);
%    size2=temp2(1);
%    
%    k=(i-1)*1/no_of_nodes;
%    
%    line([i i],[0 no_of_games],'Color','y');
%    hold on;
%    
%    for j=2:size1+1
%        plot(i,temp1(j),'r*','LineWidth',0.5);
%        hold on; 
%    end
%    
%    hold on;
%    
%    for j=2:size2+1
%        if(temp2(j)==temp1(j))
%             plot(i,temp2(j),'b*','LineWidth',0.5);
%        else
%             plot(i,temp2(j),'g*','LineWidth',0.5);
%        end
%        hold on; 
%    end    
% end