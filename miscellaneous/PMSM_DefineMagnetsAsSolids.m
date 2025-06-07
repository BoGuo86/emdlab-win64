function PMSM_DefineMagnetsAsSolids(s, mName, Nm)
for i = 1:Nm
  name = [mName, num2str(i)];
  s.defineCoil(name, 'turns', 1, 'direction', 'positive');
  s.defineWinding(name, 'extype', 'current', 'wtype', 'solid', 'exvalue', @(t) 0);
  s.addCoil2Winding(name, name);
end
end



