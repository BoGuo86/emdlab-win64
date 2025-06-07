
function [q,p] = tei2quad(t,p)

tr = triangulation(t,p);
fb = tr.freeBoundary;

e = tr.edgeAttachments(fb);

ed = tr.edges;

fb = sort(fb,2);
ed = sort(ed,2);
ed = setdiff(ed,fb,'rows')

triplot(tr)
axis off equal

end