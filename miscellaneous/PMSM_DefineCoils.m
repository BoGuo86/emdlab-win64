%% Definition of Windings
function defineWindings(s, span, PhaseA)
Ns = 3*size(PhaseA, 1);
coils = zeros(Ns,2);
for i = 1:dim_Ns
  index1 = i;
  index2 = i + span;
  if index2 > Ns
    index2 = index2 - span;
  end
  coils(i,1) = index1;
  coils(i,2) = index2;
end
