function first_blob_center(sparse_matrix)
    m, n = size(sparse_matrix)
    
    # Get row and column indices of non-zero elements
    rows, cols = findnz(sparse_matrix)[1:2]
    
    # Track visited positions to avoid counting same blob twice
    visited = Set{Tuple{Int,Int}}()
    
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
                
                return (center_row, center_col)
            end
        end
    end
    
    return nothing  # No blob found
end

function closest_min(dst_grd, tol=1e-2)
    dst_mins = sparse(dst_grd .<= tol)
    min_idc = first_blob_center(dst_mins)
    return min_idc
end

function blob_centers(sparse_matrix)
    m, n = size(sparse_matrix)
    
    # Get row and column indices of non-zero elements
    rows, cols = findnz(sparse_matrix)[1:2]
    
    # Track visited positions to avoid counting same blob twice
    visited = Set{Tuple{Int,Int}}()
    
    # Helper function for BFS considering wrap-around in column dimension
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
            
            # Check neighbors (8-connected), with wrap-around in column
            for dr in -1:1
                for dc in -1:1
                    new_row, new_col = row + dr, col + dc
                    
                    # Wrap around column if needed
                    if new_col < 1
                        new_col += n
                    elseif new_col > n
                        new_col -= n
                    end
                    
                    if 1 <= new_row <= m &&
                       sparse_matrix[new_row, new_col] == 1 &&
                       (new_row, new_col) ∉ visited
                        push!(queue, (new_row, new_col))
                    end
                end
            end
        end
        
        return blob_positions
    end
    
    # List to store centers of all blobs
    blob_centers = []
    
    # Search row by row
    for i in 1:m
        for j in 1:n
            if sparse_matrix[i, j] == 1 && (i, j) ∉ visited
                # Found start of a new blob
                blob_positions = explore_blob(i, j)
                
                # Calculate center of blob
                center_row = round(Int, mean([p[1] for p in blob_positions]))
                
                # Adjust column calculation for wrap-around
                col_positions = [p[2] for p in blob_positions]
                col_mean = mean(col_positions)
                if maximum(col_positions) - minimum(col_positions) > n / 2
                    # Adjust for wrap-around
                    col_mean = mean(mod1.(col_positions, n))
                end
                center_col = round(Int, col_mean)
                
                # Add center to the list
                push!(blob_centers, (center_row, center_col))
            end
        end
    end
    
    return blob_centers  # Return all blob centers
end

function all_mins(dst_grd, tol=1e-2)
    dst_mins = sparse(dst_grd .<= tol)
    min_idc = blob_centers(dst_mins)
    return min_idc
end