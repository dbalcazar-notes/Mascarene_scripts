rm(list=ls())

library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(GenomicRanges)
library(dplyr)

#setwd("/vast/palmer/pi/powell/DB_upload/BCF_MASC/LAPI/dsuite/try_07/sex/dsx/")
setwd("/vast/palmer/pi/powell/DB_upload/BCF_MASC/genes/rdl/")
#setwd("/vast/palmer/pi/powell/DB_upload/culex/")

ref_fasta <- "~/Reference/gen_genome/MAS/2024_06_04_MascCH02-dc_hifiasm_ph_yahs_jb_no_debris_tgs.fasta"
#ref_fasta <- "chromosome_2_extraction_reversed.fasta"
indexFa(ref_fasta)

gene <- "g3702"

bam_dir <- paste0(gene, "_reads")
out_base <- paste0(gene, "_fasta")
dir.create(out_base, showWarnings = FALSE)

bed_file <- paste0("CDS_bed_per_gene/", gene, ".bed")

bed <- read.table(
  bed_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(bed)[1:3] <- c("chr", "start0", "end")

bed <- bed %>%
  dplyr::mutate(
    start0 = as.numeric(start0),
    end = as.numeric(end),
    start = start0 + 1,
    exon_id = sprintf("%03d", dplyr::row_number()),
    len_bed = end - start0,
    len_r = end - start + 1
  )

bed

bam_files <- list.files(
  bam_dir,
  pattern = "\\.bam$",
  full.names = TRUE
)

process_bam <- function(bam_file, region_chr, region_start, region_end, exon_id) {
  
  sample_id <- sub("\\.bam$", "", basename(bam_file))
  sample_id <- sub(paste0("\\.", gene, "$"), "", sample_id)
  
  message("Procesando: ", sample_id)
  
  out_dir <- file.path(out_base, sample_id, paste0("exon_", exon_id))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  param <- ScanBamParam(
    which = GRanges(region_chr, IRanges(region_start, region_end))
  )
  
  pparam <- PileupParam(
    distinguish_nucleotides = TRUE,
    distinguish_strands = FALSE,
    include_deletions = FALSE,
    include_insertions = FALSE,
    min_base_quality = 20,
    min_mapq = 20
  )
  
  pu <- pileup(
    bam_file,
    scanBamParam = param,
    pileupParam = pparam
  )
  
  fa <- FaFile(ref_fasta)
  ref_seq <- getSeq(
    fa,
    GRanges(region_chr, IRanges(region_start, region_end))
  )
  ref_vec <- strsplit(as.character(ref_seq), "")[[1]]
  
  if (nrow(pu) == 0) {
    message("  Sin cobertura: ", sample_id)
    return(NULL)
  }
  
  ref_df <- data.frame(
    pos = region_start:region_end,
    ref = ref_vec,
    stringsAsFactors = FALSE
  )
  
  variants_all <- pu %>%
    dplyr::filter(nucleotide %in% c("A", "C", "G", "T")) %>%
    dplyr::group_by(seqnames, pos) %>%
    dplyr::mutate(depth = sum(count)) %>%
    dplyr::arrange(desc(count), .by_group = TRUE) %>%
    dplyr::summarise(
      depth = dplyr::first(depth),
      allele1 = as.character(dplyr::first(nucleotide)),
      count1 = dplyr::first(count),
      allele2 = if (dplyr::n() >= 2) as.character(nucleotide[2]) else NA_character_,
      count2 = if (dplyr::n() >= 2) count[2] else 0,
      alt_freq = count2 / depth,
      major_freq = count1 / depth,
      .groups = "drop"
    ) %>%
    dplyr::left_join(ref_df, by = "pos")
  
  fixed_alt <- variants_all %>%
    dplyr::filter(
      depth >= min_depth,
      allele1 != ref,
      major_freq >= 0.90
    )
  
  snps <- variants_all %>%
    dplyr::filter(
      depth >= min_depth,
      !is.na(allele2),
      count2 >= min_alt_count,
      alt_freq >= min_alt_freq
    )
  
  candidate_snps <- snps$pos
  
  if (length(candidate_snps) == 0) {
    
    message("  Sin SNPs heterocigotos; escribiendo consenso diploide: ", sample_id)
    
    seq_vec <- ref_vec
    
    if (nrow(fixed_alt) > 0) {
      for (i in seq_len(nrow(fixed_alt))) {
        idx <- fixed_alt$pos[i] - region_start + 1
        allele <- fixed_alt$allele1[i]
        
        if (!is.na(idx) &&
            idx >= 1 &&
            idx <= length(seq_vec) &&
            allele %in% c("A", "C", "G", "T")) {
          seq_vec[idx] <- allele
        }
      }
    }
    
    consensus_seq <- paste(seq_vec, collapse = "")
    
    dna_h1 <- DNAStringSet(consensus_seq)
    dna_h2 <- DNAStringSet(consensus_seq)
    
    names(dna_h1) <- paste0(sample_id, "_exon_", exon_id, "_hap_1")
    names(dna_h2) <- paste0(sample_id, "_exon_", exon_id, "_hap_2")
    
    writeXStringSet(
      dna_h1,
      filepath = file.path(out_dir, paste0(sample_id, "_exon_", exon_id, "_hap_1.fasta"))
    )
    
    writeXStringSet(
      dna_h2,
      filepath = file.path(out_dir, paste0(sample_id, "_exon_", exon_id, "_hap_2.fasta"))
    )
    
    write.table(
      data.frame(
        sample = sample_id,
        status = "no_heterozygous_snps",
        n_snps = 0,
        n_blocks = 0
      ),
      file = file.path(out_dir, paste0(sample_id, "_haplotypes_by_block.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    message("  OK referencia: ", sample_id)
    
    return(
      data.frame(
        sample = as.character(sample_id)[1],
        exon_id = as.character(exon_id)[1],
        chr = as.character(region_chr)[1],
        start = as.integer(region_start)[1],
        end = as.integer(region_end)[1],
        n_snps = 0,
        n_blocks = 0,
        n_fixed_alt = nrow(fixed_alt),
        out_dir = as.character(out_dir)[1],
        stringsAsFactors = FALSE
      )
    )
  }
  
  param2 <- ScanBamParam(
    which = GRanges(region_chr, IRanges(min(candidate_snps), max(candidate_snps))),
    what = c("qname", "seq", "pos", "cigar", "mapq"),
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )
  
  gal <- readGAlignments(bam_file, param = param2)
  
  if (length(gal) == 0) {
    message("  Sin reads útiles: ", sample_id)
    return(NULL)
  }
  
  seqs <- mcols(gal)$seq
  read_names <- mcols(gal)$qname
  names(seqs) <- read_names
  
  read_alleles <- matrix(
    NA,
    nrow = length(gal),
    ncol = length(candidate_snps),
    dimnames = list(read_names, as.character(candidate_snps))
  )
  
  get_base_from_read <- function(gal_i, seq_i, snp_pos) {
    
    ref_ranges <- cigarRangesAlongReferenceSpace(
      cigar(gal_i),
      pos = start(gal_i),
      ops = c("M", "=", "X")
    )[[1]]
    
    query_ranges <- cigarRangesAlongQuerySpace(
      cigar(gal_i),
      ops = c("M", "=", "X")
    )[[1]]
    
    if (length(ref_ranges) == 0) return(NA)
    
    hits <- which(snp_pos >= start(ref_ranges) & snp_pos <= end(ref_ranges))
    if (length(hits) == 0) return(NA)
    
    j <- hits[1]
    offset <- snp_pos - start(ref_ranges)[j]
    q_pos <- start(query_ranges)[j] + offset
    
    base <- substr(as.character(seq_i), q_pos, q_pos)
    
    if (base %in% c("A", "C", "G", "T")) base else NA
  }
  
  read_alleles <- matrix(
    NA,
    nrow = length(gal),
    ncol = length(candidate_snps),
    dimnames = list(names(gal), as.character(candidate_snps))
  )
  
  for (i in seq_along(gal)) {
    for (snp in candidate_snps) {
      read_alleles[i, as.character(snp)] <- get_base_from_read(
        gal[i],
        seqs[i],
        snp
      )
    }
  }
  
  # Esta matriz conserva todos los reads, incluso si cubren 1 solo SNP.
  read_alleles_all <- read_alleles
  
  # Esta matriz se usa solo para definir fase/conexión entre SNPs.
  read_alleles_phase <- read_alleles[
    rowSums(!is.na(read_alleles)) >= 2,
    ,
    drop = FALSE
  ]
  
  if (nrow(read_alleles_phase) == 0) {
    
    message("  SNPs no conectados por reads; usando SNPs individuales no faseados: ", sample_id)
    
    haplotypes <- snps %>%
      arrange(pos) %>%
      mutate(
        block = seq_len(n()),
        position = pos,
        haplotype_1 = as.character(allele1),
        haplotype_2 = as.character(allele2)
      ) %>%
      dplyr::select(
        block,
        position,
        haplotype_1,
        haplotype_2
      )
    
  } else {
    
    snps_chr <- colnames(read_alleles_phase)
    
    adj <- setNames(
      vector("list", length(snps_chr)),
      snps_chr
    )
    
    for (i in seq_len(nrow(read_alleles_phase))) {
      
      covered <- names(which(!is.na(read_alleles_phase[i, ])))
      
      if (length(covered) >= 2) {
        
        combs <- combn(covered, 2)
        
        for (k in seq_len(ncol(combs))) {
          
          a <- combs[1, k]
          b <- combs[2, k]
          
          adj[[a]] <- unique(c(adj[[a]], b))
          adj[[b]] <- unique(c(adj[[b]], a))
        }
      }
    }
    
    visited <- setNames(
      rep(FALSE, length(snps_chr)),
      snps_chr
    )
    
    blocks <- list()
    
    for (s in snps_chr) {
      
      if (!visited[s]) {
        
        queue <- s
        block <- c()
        
        while (length(queue) > 0) {
          
          x <- queue[1]
          queue <- queue[-1]
          
          if (!visited[x]) {
            
            visited[x] <- TRUE
            
            block <- c(block, x)
            
            neigh <- adj[[x]]
            
            queue <- c(
              queue,
              neigh[!visited[neigh]]
            )
          }
        }
        
        blocks[[length(blocks) + 1]] <- sort(
          as.numeric(block)
        )
      }
    }
    
    snps_in_blocks <- sort(
      unique(unlist(blocks))
    )
    
    isolated_snps <- setdiff(
      candidate_snps,
      snps_in_blocks
    )
    
    if(length(isolated_snps) > 0){
      
      for(s in isolated_snps){
        blocks[[length(blocks)+1]] <- s
      }
    }
    
    blocks <- blocks[
      order(sapply(blocks, min))
    ]
    
    phase_block_by_linkage <- function(mat_use, min_link_reads = 1) {
      
      snp_names <- colnames(mat_use)
      n_snps <- ncol(mat_use)
      
      hap1 <- rep(NA_character_, n_snps)
      hap2 <- rep(NA_character_, n_snps)
      names(hap1) <- snp_names
      names(hap2) <- snp_names
      
      # SNP inicial: solo para anclar el bloque
      first_snp <- snp_names[1]
      
      alleles_first <- mat_use[, first_snp]
      alleles_first <- alleles_first[alleles_first %in% c("A", "C", "G", "T")]
      tab_first <- sort(table(alleles_first), decreasing = TRUE)
      
      if (length(tab_first) >= 2) {
        hap1[first_snp] <- names(tab_first)[1]
        hap2[first_snp] <- names(tab_first)[2]
      } else if (length(tab_first) == 1) {
        hap1[first_snp] <- names(tab_first)[1]
        hap2[first_snp] <- names(tab_first)[1]
      }
      
      # Propagar fase solo si hay reads conectando SNPs
      for (snp in snp_names[-1]) {
        
        phased_snps <- snp_names[!is.na(hap1) & !is.na(hap2)]
        
        support_h1 <- character()
        support_h2 <- character()
        
        for (prev_snp in phased_snps) {
          
          reads_link <- !is.na(mat_use[, prev_snp]) & !is.na(mat_use[, snp])
          
          if (sum(reads_link) < min_link_reads) next
          
          df_link <- data.frame(
            prev = mat_use[reads_link, prev_snp],
            curr = mat_use[reads_link, snp],
            stringsAsFactors = FALSE
          ) %>%
            dplyr::filter(
              prev %in% c("A", "C", "G", "T"),
              curr %in% c("A", "C", "G", "T")
            )
          
          support_h1 <- c(
            support_h1,
            df_link$curr[df_link$prev == hap1[prev_snp]]
          )
          
          support_h2 <- c(
            support_h2,
            df_link$curr[df_link$prev == hap2[prev_snp]]
          )
        }
        
        h1_tab <- sort(table(support_h1), decreasing = TRUE)
        h2_tab <- sort(table(support_h2), decreasing = TRUE)
        
        if (length(h1_tab) > 0) {
          hap1[snp] <- names(h1_tab)[1]
        }
        
        if (length(h2_tab) > 0) {
          hap2[snp] <- names(h2_tab)[1]
        }
      }
      
      hap1[is.na(hap1)] <- "N"
      hap2[is.na(hap2)] <- "N"
      
      list(
        haplotype_1 = hap1,
        haplotype_2 = hap2
      )
    }
    
    hap_out <- list()
    
    for (b in seq_along(blocks)) {
      
      block_snps <- as.character(blocks[[b]])
      
      mat <- read_alleles_all[
        ,
        block_snps,
        drop = FALSE
      ]
      
      mat_use <- mat[
        rowSums(!is.na(mat)) >= 1,
        ,
        drop = FALSE
      ]
      
      phased <- phase_block_by_linkage(mat_use)
      
      pat1 <- phased$haplotype_1
      pat2 <- phased$haplotype_2
      
      hap_out[[b]] <- data.frame(
        block = b,
        position = as.integer(block_snps),
        haplotype_1 = pat1,
        haplotype_2 = pat2,
        stringsAsFactors = FALSE
      )
    }
    
    haplotypes <- do.call(rbind, hap_out)
  }
  
  write.table(
    haplotypes,
    file = file.path(
      out_dir,
      paste0(sample_id, "_haplotypes_by_block.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  make_sequence_from_block <- function(block_id, hap_col) {
    
    seq_vec <- ref_vec
    
    ## Primero aplicar cambios fijados a ambos haplotipos
    if (nrow(fixed_alt) > 0) {
      for (i in seq_len(nrow(fixed_alt))) {
        idx <- fixed_alt$pos[i] - region_start + 1
        allele <- fixed_alt$allele1[i]
        
        if (!is.na(idx) &&
            idx >= 1 &&
            idx <= length(seq_vec) &&
            allele %in% c("A", "C", "G", "T")) {
          seq_vec[idx] <- allele
        }
      }
    }
    
    ## SNPs heterocigotos del bloque actual
    block_df <- haplotypes %>%
      dplyr::filter(block == block_id)
    
    ## SNPs heterocigotos de otros bloques quedan como N
    other_snp_pos <- setdiff(
      haplotypes$position,
      block_df$position
    )
    
    if (length(other_snp_pos) > 0) {
      other_idx <- other_snp_pos - region_start + 1
      other_idx <- other_idx[other_idx >= 1 & other_idx <= length(seq_vec)]
      seq_vec[other_idx] <- "N"
    }
    
    ## SNPs del bloque actual reciben el alelo haplotípico
    for (i in seq_len(nrow(block_df))) {
      idx <- block_df$position[i] - region_start + 1
      allele <- block_df[[hap_col]][i]
      
      if (!is.na(idx) &&
          idx >= 1 &&
          idx <= length(seq_vec) &&
          allele %in% c("A", "C", "G", "T")) {
        seq_vec[idx] <- allele
      } else if (!is.na(idx) &&
                 idx >= 1 &&
                 idx <= length(seq_vec)) {
        seq_vec[idx] <- "N"
      }
    }
    
    paste(seq_vec, collapse = "")
  }
  
  for (b in unique(haplotypes$block)) {
    
    seq_h1 <- make_sequence_from_block(b, "haplotype_1")
    seq_h2 <- make_sequence_from_block(b, "haplotype_2")
    
    dna_h1 <- DNAStringSet(seq_h1)
    dna_h2 <- DNAStringSet(seq_h2)
    
    names(dna_h1) <- paste0(sample_id, "_exon_", exon_id, "_BLOCK_", b, "_hap_1")
    names(dna_h2) <- paste0(sample_id, "_exon_", exon_id, "_BLOCK_", b, "_hap_2")
    
    writeXStringSet(
      dna_h1,
      filepath=file.path(
        out_dir,
        paste0(
          sample_id,
          "_exon_",
          exon_id,
          "_BLOCK_",
          b,
          "_hap_1.fasta"
        )
      )
    )
    
    writeXStringSet(
      dna_h2,
      filepath=file.path(
        out_dir,
        paste0(
          sample_id,
          "_exon_",
          exon_id,
          "_BLOCK_",
          b,
          "_hap_2.fasta"
        )
      )
    )
  }
  
  message("  OK: ", sample_id)
  
  return(
    data.frame(
      sample = as.character(sample_id)[1],
      exon_id = as.character(exon_id)[1],
      chr = as.character(region_chr)[1],
      start = as.integer(region_start)[1],
      end = as.integer(region_end)[1],
      n_snps = nrow(haplotypes),
      n_blocks = length(unique(haplotypes$block)),
      out_dir = as.character(out_dir)[1],
      stringsAsFactors = FALSE
    )
  )
}

min_depth <- 8
min_alt_count <- 3
min_alt_freq <- 0.25

summary_list <- list()

counter <- 1

for (e in seq_len(nrow(bed))) {
  
  region_chr <- bed$chr[e]
  region_start <- bed$start[e]
  region_end <- bed$end[e]
  exon_id <- bed$exon_id[e]
  
  message("========== EXON ", exon_id, ": ",
          region_chr, ":", region_start, "-", region_end,
          " ==========")
  
  for (bam in bam_files) {
    
    summary_list[[counter]] <- process_bam(
      bam_file = bam,
      region_chr = region_chr,
      region_start = region_start,
      region_end = region_end,
      exon_id = exon_id
    )
    
    counter <- counter + 1
  }
}

bed
