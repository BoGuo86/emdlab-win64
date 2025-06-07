classdef ali_class
  properties
    a
  end
  methods
    function obj = ali_class()
      a=MException('','ag');
      throw(a)
    end
  end
end
