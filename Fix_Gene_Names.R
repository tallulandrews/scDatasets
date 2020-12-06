fix_gene_names <- function(curr_names, ref_names) {
	curr_names <- strsplit(unlist(curr_names), "[;,:]")
	best_names <- lapply(curr_names, function(x) {
		if (sum(grepl("Mar-", x))>0) {
			x <- c(x, sub("Mar-", "MARCH", x), sub("Mar-", "March", x))
		}
		if (sum(grepl("-Mar", x))>0) {
			x <- c(x, sub("([1234567890]+)-Mar", "MARCH\\1", x), sub("([1234567890]+)-Mar", "March\\1", x))
		}
		if (sum(grepl("Dec-", x))>0) {
			x <- c(x, sub("Dec-", "DEC", x), sub("Dec-", "Dec", x))
		}
		if (sum(grepl("-Dec", x))>0) {
			x <- c(x, sub("([1234567890]+)-Dec", "DEC\\1", x), sub("([1234567890]+)-Dec", "Dec\\1", x))
		}
		if (sum(grepl("Sep-", x))>0) {
			x <- c(x, sub("Sep-", "SEPT", x), sub("Sep-", "Sept", x))
		}
		if (sum(grepl("-Sep", x))>0) {
			x <- c(x, sub("([1234567890]+)-Sep", "SEPT\\1", x), sub("([1234567890]+)-Sep", "Sept\\1", x))
		}
		

		if (sum(x %in% ref_names) >=1) {
			x <- x[x %in% ref_names]
		}
		return(x[1]);
	})
	return(unlist(best_names))
}
