
function t = tds_getTree()

% root node
root = javax.swing.tree.DefaultMutableTreeNode();

% main
superNode = javax.swing.tree.DefaultMutableTreeNode("Main");
root.add(superNode);
nodes = ["Form", "N parts"];
for i = 1:length(nodes)
    superNode.add(javax.swing.tree.DefaultMutableTreeNode(nodes(i)));
end

% geometry
superNode = javax.swing.tree.DefaultMutableTreeNode("Geometry");
root.add(superNode);
nodes = ["Form", "DS", "HF", "ES", "ESR", "BS", "SS", "BJ"];
for i = 1:length(nodes)
    superNode.add(javax.swing.tree.DefaultMutableTreeNode(nodes(i)));
end

% physics
superNode = javax.swing.tree.DefaultMutableTreeNode("Physics");
root.add(superNode);
nodes = ["Sigma", "E-Modul"];
for i = 1:length(nodes)
    superNode.add(javax.swing.tree.DefaultMutableTreeNode(nodes(i)));
end

% windings
superNode = javax.swing.tree.DefaultMutableTreeNode("Windings");
root.add(superNode);
nodes = ["uk_ind", "uk-ohm", "N"];
for i = 1:length(nodes)
    superNode.add(javax.swing.tree.DefaultMutableTreeNode(nodes(i)));
end

% winding parts
superNode = javax.swing.tree.DefaultMutableTreeNode("Parts");
root.add(superNode);
nodes = ["part1", "part2"];
for i = 1:length(nodes)
    superNode.add(javax.swing.tree.DefaultMutableTreeNode(nodes(i)));
end

% ms2df
superNode = javax.swing.tree.DefaultMutableTreeNode("MS2DF");
root.add(superNode);
nodes = ["Nh", "part2"];
for i = 1:length(nodes)
    superNode.add(javax.swing.tree.DefaultMutableTreeNode(nodes(i)));
end

t = javax.swing.JTree(root);
t.setRootVisible(false);
t.setShowsRootHandles(true);
t.getSelectionModel().setSelectionMode(javax.swing.tree.TreeSelectionModel.SINGLE_TREE_SELECTION);

end