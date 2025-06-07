classdef emdlab_mcs_circularShell < handle & g2d_constants
  properties (SetAccess = private)
    % inner points
    ips (:,2) double;
    % outer points
    ops (:,2) double;
    % inner point angles
    ipas (:,1) double;
    % outer point angles
    opas (:,1) double;
    % index of inner in outer
    ii2o (:,2) double;
    % index of outer in inner
    io2i (:,2) double;
    % shell center
    csh (1,2) double;
    % shell radius
    rsh (1,1) double;
    % index of sorted points
    iosp (:,1) double;
  end
  properties (Dependent = true)
    % number of inner points
    Nps (1,1) double;
  end
  methods
    function obj = emdlab_mcs_circularShell(ps)
      % evaluation of shell center
      obj.csh = sum(ps)/size(ps,1);
      % assigning points to inner and outer
      obj.ips = ps;
      obj.ops = ps; 
      % evaluation of shell radius
      tmp = [ps(:,1)-obj.csh(1), ps(:,2)-obj.csh(2)];
      tmp = sqrt(sum(tmp.^2, 2));
      obj.rsh = mean(tmp);
      if sum(abs(tmp-obj.rsh)) > obj.Nps*obj.gleps
        error('points do not form a cirlcular shell.');
      end      
      obj.updateData;
    end
    function y = get.Nps(obj)
      y = size(obj.ips, 1);
    end
    function updateData(obj)
      % inner points
      obj.ipas = atan_02pi([obj.ips(:,1) - obj.icc(1), obj.ips(:,2) - obj.icc(2)]);
      % outer points
      obj.opas = atan_02pi([obj.ops(:,1) - obj.occ(1), obj.ops(:,2) - obj.occ(2)]);
    end
    function obj = rotateInner(obj, angle)
      obj.ips = ext_protate2(obj.ips, angle);
      obj.updateData;
    end
    function obj = rotateOuter(obj, angle)
      obj.ops = ext_protate2(obj.ops, angle);
      obj.updateData;
    end
  end
end
