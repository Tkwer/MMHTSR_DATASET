% This function grows clusters by using BFS algorithm
function [center, result_max] = grow_clusters(pointarray)
    result_max = []; % initialize the variable to store the largest cluster
    while ~isempty(pointarray) % if the pointarray is not empty
        result = BFS(pointarray,pointarray(1,:)); % apply BFS algorithm to get one cluster
        size_w = size(result);
        size_w_max = size(result_max);
        if size_w(1)>=size_w_max(1) % if the current cluster is larger than the largest cluster found so far
            result_max = result; % update the largest cluster
        end
        % remove the points in the current cluster from the pointarray
        for index =1:size_w(1)
            rowsToRemove = ismember(pointarray, result(index, :), 'rows');
            pointarray(rowsToRemove, :) = [];
        end
    end
    % calculate the center of the largest cluster
    center = sum(result_max,1)./length(result_max(:,1));
end

% This function implements the BFS algorithm to search for a cluster
function result = BFS(pointarray,s)
    queue = [];  % create a queue
    queue = [queue; s]; % add the starting point to the queue
    seen = [];  % record the visited points
    seen = [seen; s];
    while ~isempty(queue) % while the queue is not empty
        vertex = queue(1,:); % get the first point in the queue
        queue(1,:)=[];% dequeue
        k = get8Neighborhood(vertex); % get the 8 neighbors of the current point
        index_ = ismember(k,pointarray,'rows');
        kk =[];
        % remove the points that are not in the pointarray
        for j=1:9
            if index_(j)==1&&j~=5 % the current point is always in the pointarray, so we don't need to update it
                kk = [kk;k(j,:)];
            end
        end
        next = kk; % get the neighbors that are in the pointarray
        size_2 = size(next);
        for i =1:size_2(1) % for each neighbor
            if ~ismember(next(i,:),seen,'rows') % if it has not been visited
                queue = [queue; next(i,:)]; % enqueue
                seen = [seen;next(i,:)]; % mark as visited
            end
        end
    end
    result = seen; % return the visited points as a cluster
end

% This function returns the 8 neighbors of a center point
function point=get8Neighborhood(centerpoint)
    point = zeros(9,2); % initialize the neighbor points
    for i=1:3
        for j=1:3
             point((i-1)*3+j,1) = centerpoint(1)+i-2; % get the x coordinate of the neighbor point
             point((i-1)*3+j,2) = centerpoint(2)+j-2; % get the y coordinate of the neighbor point
        end
    end
end
