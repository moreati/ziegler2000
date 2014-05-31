function [ConsecRoots,t] = FindConsecRoots(BaseRoots)

if length(BaseRoots)<2
   ConsecRoots = -1;
   t = -1;
   return
end

derivative = diff(BaseRoots);         % first derivative
CRI = find(derivative==1);
ConsecRoots = BaseRoots(CRI);
DiffCount = 1;

% Take higher derivatives to separate clusters of consecutve roots
% and find the largest group
while length(CRI)>1
   DiffCount = DiffCount+1;
   if DiffCount>20
      disp 'Error: too many derivatives!'
      dbcont
   end
   derivative = diff(ConsecRoots);
   CRI = find(derivative==1);
   ConsecRoots = ConsecRoots(CRI);
end

ConsecRoots = [ConsecRoots:ConsecRoots+DiffCount];
t = length(ConsecRoots)/2;




