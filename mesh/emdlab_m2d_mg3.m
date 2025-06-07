%% EMDLAB
% mesh generator:
% this mesh generator generates mesh with two constraints:
% 1) maximum allowaible length for the edges
% 2) aspect ratio of triangles

function mz = emdlab_m2d_mg3(f, v, f1, v1)

    emdlab_g2d_validateFV(f, v);

    % evaluation of max length according to boundary edges
    tmp = v(f(:, 2), :) - v(f(:, 1), :);
    tmp = sqrt(sum(tmp.^2, 2));
    maxLength = max(tmp);

    % setting aspect ratio
    aspectRatio = 2;

    % length of boundary edges
    el = sqrt(sum((v(f(:,1),:)-v(f(:,2),:)).^2,2));

    % mid point of edges
    midp = (v(f(:,1),:)+v(f(:,2),:))/2;

    % edge length function
    fh = scatteredInterpolant(midp(:,1),midp(:,2),el);
    fh = @(p) fh(p(:,1),p(:,2));

    % first delaunay triangulation
    dt = delaunayTriangulation(v, f);
    mz = emdlab_m2d_tmz(dt.ConnectivityList(ext_dpoly2d(dt.incenter, f1, v1) <- 1e-6, :), dt.Points);

    % number of iterations
    iter = 0;

    while iter < 50

        ff = mz.edges(mz.bedges, :);
        tmp = mz.getCenterOfElements;
%         tmp = tmp(mz.getMaxEdgeLength > maxLength | mz.getAspectRatio > aspectRatio, :);
        tmp = tmp(mz.getMaxEdgeLength > 1.2*fh(tmp) & mz.getAspectRatio > aspectRatio, :);

        if ~any(tmp)
            break
        end

        vv = [mz.nodes; tmp];
        dt = delaunayTriangulation(vv, ff);
        mz = emdlab_m2d_tmz(dt.ConnectivityList(ext_dpoly2d(dt.incenter, f1, v1) <- 1e-6, :), dt.Points);
        mz.moveNodes;

        % next iteration
        iter = iter + 1;

    end

end
