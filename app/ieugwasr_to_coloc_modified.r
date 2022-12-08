ieugwasr_to_coloc_custom <- function(id1, eqtl_df, chrompos, type1=NULL, type2=NULL)
{
    id2 <- 'eQTL-Gen'
	tab1 <- ieugwasr::associations(id=id1, variants=chrompos) %>% subset(., !duplicated(rsid))
	# tab2 <- ieugwasr::associations(id=id2, variants=chrompos) %>% subset(., !duplicated(rsid))
	tab2 <- eqtl_df %>% subset(., !duplicated(rsid)) # tab2 is the eQTL dataset
	commonsnps <- tab1$rsid[tab1$rsid %in% tab2$rsid]
	tab1 <- tab1[tab1$rsid %in% commonsnps, ] %>% dplyr::arrange(rsid)
	tab2 <- tab2[tab2$rsid %in% commonsnps, ] %>% dplyr::arrange(rsid)
	stopifnot(all(tab1$rsid == tab2$rsid))

	index <- as.character(tab1$ea) == as.character(tab2$ea) &
			as.character(tab1$nea) == as.character(tab2$nea) &
			as.character(tab1$rsid) == as.character(tab2$rsid) &
			tab1$position == tab2$position
	stopifnot(sum(index) > 0)
	tab1$eaf <- as.numeric(tab1$eaf)
	tab2$eaf <- as.numeric(tab2$eaf)
	tab1$eaf[which(tab1$eaf > 0.5)] <- 1 - tab1$eaf[which(tab1$eaf > 0.5)]
	tab2$eaf[which(tab2$eaf > 0.5)] <- 1 - tab2$eaf[which(tab2$eaf > 0.5)]
	s <- sum(is.na(tab1$eaf))
	if(s > 0)
	{
		warning(s, " out of ", nrow(tab1), " variants have missing allele frequencies in ", id1, ". Setting to 0.5")
		tab1$eaf[is.na(tab1$eaf)] <- 0.5
	}
	s <- sum(is.na(tab2$eaf))
	if(s > 0)
	{
		warning(s, " out of ", nrow(tab2), " variants have missing allele frequencies in ", id2, ". Setting to 0.5")
		tab2$eaf[is.na(tab2$eaf)] <- 0.5
	}

	info1 <- ieugwasr::gwasinfo(id1)
	type1 <- get_type(info1, type1)
	# info2 <- ieugwasr::gwasinfo(id2)
	type2 <- "quant"


	tab1$n[is.na(tab1$n)] <- info1$sample_size
    # swap any NA values in eqtl table w/ median sample size
    tab2$n[is.na(tab2$n)] <- median(tab2$n)

	tab1 <- tab1[index,] %>% {list(pvalues = .$p, N = .$n, MAF = .$eaf, beta = .$beta, varbeta = .$se^2, type = type1, snp = .$rsid, z = .$beta / .$se, chr = .$chr, pos = .$position, id = id1)}
	tab2 <- tab2[index,] %>% {list(pvalues = .$p, N = .$n, MAF = .$eaf, 
        # beta = NULL, varbeta = NULL, 
        type = type2, snp = .$rsid, z = .$z, chr = .$chr, pos = .$position, id = id2)}

	if(type1 == "cc")
	{
		tab1$s <- info1$ncase / info1$sample_size
	}

	if(type2 == "cc")
	{
		tab2$s <- info2$ncase / info2$sample_size
	}

	return(list(dataset1=tab1, dataset2=tab2))
}

get_type <- function(info, typex)
{
	if(!is.null(typex))
	{
		stopifnot(typex %in% c("cc", "quant"))
		return(typex)
	} else if(is.na(info$unit)) {
		if(! "ncase" %in% names(info))
		{
			info$ncase <- NA
		}
		if(is.na(info$ncase))
		{
			message("Type information not available for ", info$id, ". Assuming 'quant' but override using 'type' arguments.")
			return("quant")			
		} else {
			message("No units available but assuming cc due to number of cases being stated")
			return("cc")
		}
	} else {
		return(ifelse(info$unit %in% c("logOR", "log odds"), "cc", "quant"))
	}
}