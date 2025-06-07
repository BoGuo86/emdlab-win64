
f = figure('menu', 'none');

[jp, wp] = javacomponent(javax.swing.JPanel, [], f);
set(wp, 'units', 'norm', 'position', [0,0,1,1]);


jp.setLayout(java.awt.BorderLayout)

p = javax.swing.JPanel;
p.setPreferredSize(java.awt.Dimension(10,100))
jp.add(p, java.awt.BorderLayout.NORTH)


jp.setLayout(java.awt.FlowLayout)

t = javax.swing.JToolBar;
t.setPreferredSize(java.awt.Dimension(10,25))
t1 = javax.swing.JToolBar;
t1.setPreferredSize(java.awt.Dimension(10,25))

p.add(t1)
p.add(t)