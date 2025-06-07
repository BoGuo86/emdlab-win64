
function test_java

f = figure('HandleVisibility', 'off');
b1 = javax.swing.JButton('click me');
t1 = javax.swing.JTextField('5');
t2 = javax.swing.JTextField('8');
t3 = javax.swing.JTextField('8');
javacomponent(b1, [50,50,80,20], f);
 javacomponent(t1, [150,50,80,20], f);
  javacomponent(t2, [150,100,80,20], f);
    javacomponent(t3, [150,150,80,20], f);
% Set the Matlab callbacks to the JButton's events
set(b1, 'MouseEnteredCallback', @(h,e)b1.setForeground(java.awt.Color.RED))
set(b1, 'MouseExitedCallback',  @(h,e)b1.setForeground(java.awt.Color.BLUE))
set(b1, 'MouseClickedCallback',  @myf)


function myf(h,e)

a = str2double(t1.getText) +str2double(t2.getText);
t3.setText(num2str(a));
end


end