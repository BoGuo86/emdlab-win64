function y = getPQUANTITY(qName)
switch qName
    case 'length'
        y = PQUANTITYC('m');
        y.addUnit('um', 1e-6);
        y.addUnit('mm', 1e-3);
        y.addUnit('km', 1e3);
    case 'power'
        y = PQUANTITYC('w');
        y.addUnit('mw', 1e-3);
        y.addUnit('kw', 1e3);
    case 'magnticFlux'
        y = PQUANTITYC('Tesla');
    case 'angularSpeed'
        y = PQUANTITYC('rad/s');
        y.addUnit('rpm', pi/30);
    case 'none'
        y = PQUANTITYC('none');
    otherwise
        error('Invalid Quantity type');
end
end
