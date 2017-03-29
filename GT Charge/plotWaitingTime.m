clear;

no_of_games=textread('Put_NO_OF_GAMES.txt','%d');
% disp(no_of_games);

no_of_nodes=textread('Put_NO_OF_NODES.txt','%d');
% disp(no_of_nodes);

axis([0 no_of_nodes 0 no_of_games]);
hold on;


waitingTime=textread('Put_AVERAGE_WAITING.txt','%f');
deadTime=textread('Put_DEAD_TIME.txt','%d');
disp(waitingTime(1));

waitingTimeModified=textread('F:\\project8\\GT Charge Modified\\Put_AVERAGE_WAITING.txt','%f');
deadTimeModified=textread('F:\\project8\\GT Charge Modified\\Put_DEAD_TIME.txt','%d');

%GTCharge
% plot(waitingTime,'g*','LineWidth',0.5);
% hold on;

% for i=1:no_of_nodes
% %     if deadTime(i)<= no_of_games
% %         line([i i],[0 no_of_games],'Color','y');
% %     end
% plot(deadTime,'b*','LineWidth',0.5);
% end
% hold on;


%GTCharge Modified
% plot(waitingTimeModified,'b*','LineWidth',0.5);
% hold on;

for i=1:no_of_nodes
%     if deadTimeModified(i)<= no_of_games
%         line([i i],[0 no_of_games],'Color','y');
%     end
end
plot(deadTimeModified,'r*','LineWidth',0.5);
hold on;