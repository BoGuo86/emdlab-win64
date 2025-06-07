% find whether given graph is star or not
function testGraph()
clc;
% example 1
g = [0,1;1,0];
isStarGraph(g)
% example 2
g = [0,1,0;1,0,0;1,1,1];
isStarGraph(g)
% example 3
g = [0,1,1,1;1,0,0,0;1,0,0,0;1,0,0,0];
isStarGraph(g)
% example 4
g = [0,0,1,0;1,0,0,0;1,0,0,0;1,0,0,0];
isStarGraph(g)
% example 5
g = [0,0,1,0;0,0,1,0;1,1,0,1;0,0,1,0];
isStarGraph(g)
end
% print function
function isStarGraph(g)
disp('**************************************************');
disp('your graph is:');
disp(g);
% displaying the result
if checkStar(g)
  disp('your graph is a star graph.');
else
  disp('your graph is not a star graph.');
end
disp('**************************************************');
end
% checkStar function
function y = checkStar(g)
% getting size of graph nodes
[m,n] = size(g);
% check for irregular graphs
if m ~= n
  error('graph must be an squer matrix.');
end
% check for single node graphs
if m == 1
  if g == 0
    y = true; return;
  else
    y = false; return;
  end
  % check for multiple node graphs
elseif m>1
  % loop over nodes
  for i = 1:m
    % testing ith node
    if g(i,i) == 0
      if all(g(i,1:i-1)) && all(g(i,i+1:end)) && all(g(1:i-1,i)) && all(g(i+1:end,i)) ...
          && all(all(~g(1:i-1,1:i-1))) && all(all(~g(i+1:end,i+1:end))) ...
          && all(all(~g(i+1:end,1:i-1))) && all(all(~g(1:i-1,i+1:end)))
        y = true; return;
      else
        % going to begginig of loop for testing another node
        if i<m
          continue;
        else
          y = false; return;
        end
      end
    else
      y = false; return;
    end
  end
else
  % check for empty graphs
  error('graph is empty.');
end
end
