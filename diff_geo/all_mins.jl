function all_blob_centers(sparse_matrix)
    m, n = size(sparse_matrix)
    
    # Get row and column indices of non-zero elements
    rows, cols = findnz(sparse_matrix)[1:2]
    
    # Track visited positions to avoid counting same blob twice
    visited = Set{Tuple{Int,Int}}()
    centers = []  # List to store centers of all blobs
    
    # Helper function for BFS
    function explore_blob(start_row, start_col)
        blob_positions = Tuple{Int,Int}[]
        queue = [(start_row, start_col)]
        
        while !isempty(queue)
            row, col = popfirst!(queue)
            
            if (row, col) in visited
                continue
            end
            
            push!(visited, (row, col))
            push!(blob_positions, (row, col))
            
            # Check neighbors (8-connected)
            for dr in -1:1
                for dc in -1:1
                    new_row, new_col = row + dr, col + dc
                    if 1 <= new_row <= m && 1 <= new_col <= n &&
                       sparse_matrix[new_row, new_col] == 1 &&
                       (new_row, new_col) ∉ visited
                        push!(queue, (new_row, new_col))
                    end
                end
            end
        end
        
        return blob_positions
    end
    
    # Search row by row
    for i in 1:m
        for j in 1:n
            if sparse_matrix[i, j] == 1 && (i, j) ∉ visited
                # Found start of a new blob
                blob_positions = explore_blob(i, j)
                
                # Calculate center of blob
                center_row = round(Int, mean([p[1] for p in blob_positions]))
                center_col = round(Int, mean([p[2] for p in blob_positions]))
                
                push!(centers, (center_row, center_col))
            end
        end
    end
    
    return centers  # Return list of all blob centers
end

function all_mins(dst_grd, tol=1e-2)
    dst_mins = sparse(dst_grd .<= tol)
    min_centers = all_blob_centers(dst_mins)
    return min_centers
end