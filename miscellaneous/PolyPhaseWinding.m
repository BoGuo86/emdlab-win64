classdef PolyPhaseWinding
  properties
    Nphases (1,1) double {mustBePositive, mustBeInteger} = 1;
    Nslots (1,1) double {mustBePositive, mustBeInteger} = 1;
    Npoles (1,1) double {mustBePositive, mustBeInteger} = 1;
    Nlyrs (1,1) double {mustBePositive, mustBeInteger} = 1;
    Nfraq (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
  end
  properties (Dependent = true)
    q
  end
  methods
    function obj = PolyPhaseWinding(Nphases, Nslots, Npoles, Nlyrs, Nfraq)
      obj.Nphases = Nphases;
      obj.Nslots = Nslots;
      obj.Npoles = Npoles;
      obj.Nlyrs = Nlyrs;
      obj.Nfraq = Nfraq;
      % check interrelationships
      if rem(Nslots, Nphases)
        error('number of slots must be a multiple number of phases.');
      end
      if rem(Npoles, 2)
        error('number of poles must be an even number.');
      end
      if ~ismember(Nlyrs, [1,2])
        error('number of layers must be 1 or 2.');
      end
      
    end
    function y = get.q(obj)
      y = obj.Ns/obj.p/obj.m;
    end
    function y = getPhases(obj)
      if obj.Lyrs == 1
        for i = 1:obj.m
          y.(['Phase',num2str(i)]).pos = [];
          y.(['Phase',num2str(i)]).neg = [];
          for j = 1:obj.p/2
            for k = 1:obj.q
              y.(['Phase',num2str(i)]).pos(end+1) = (j-1)*2*obj.m*obj.q+k;
              y.(['Phase',num2str(i)]).neg(end+1) = (2*j-1)*obj.m*obj.q+k;
            end
          end
        end
      end
    end
  end
end
