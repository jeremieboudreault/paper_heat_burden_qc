# create_lagged_mat.R

create_lagged_mat <- function(var) {
    
    # Loop on all RSS.
    mat_l <- lapply(rss_sel, function(rss_i) {
        
        return(as.matrix(merged_data[
            i    = (RSS == rss_i & MONTH %in% months),
            j    = paste0(var, 0:max_lag), 
            with = FALSE
        ]))
        
    })
    
    # Add names.
    names(mat_l) <- rss_sel
    
    # Return.
    return(mat_l)
    
}
