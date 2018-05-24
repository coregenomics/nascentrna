genes_overlapping <- GenomicRanges::GRanges(
    c(
        "chr1:10001-15000:+",
        "chr1:15000-20000:+",       # Overlaps with gene above.
        "chr1:30001-40000:-"
    ),
    gene_id = 1:3)
genes <- GenomicRanges::GRanges(
    c(
        "chr1:10001-15000:+",
        "chr1:30001-40000:-"
    ),
    gene_id = 1:2)
tss <- GenomicRanges::GRanges(
    c(
        ## TSSs related to gene 1 chr:10000-15000+ (# nolint)
        "chr1:9801-9900:-",             # uaRNA.
        "chr1:9901-10000:-",            # uaRNA overlaps with promoter.
        "chr1:9900-10100:-",            # aRNA overlaps with promoter.
        "chr1:10001-10200:-",           # daTSS that overlaps with promoter.
        "chr1:10101-10300:+",           # Alternative TSS at 100 bp.
        "chr1:10101-10300:-",           # daTSS at 100 bp.
        "chr1:10301-10400:-",           # daTSS at 300 bp.
        "chr1:10501-10600:-",           # daRNA at 500 bp.
        ## TSSs related to gene 2 chr1:30000-40000- (# nolint)
        "chr1:29800-30000:-",           # 3' convergent RNA.
        "chr1:40101-40200:+",           # uaRNA.
        "chr1:40001-40100:+",           # uaRNA overlaps with promoter.
        "chr1:39901-40100:+",           # aRNA overlaps with promoter.
        "chr1:39801-40000:+",           # daTSS that overlaps with promoter.
        "chr1:39701-39900:-",           # Alternative TSS at 100 bp.
        "chr1:39701-39900:+",           # daTSS at 100 bp.
        "chr1:39601-39700:+",           # daTSS at 300 bp.
        "chr1:39401-39500:+"            # daRNA at 500 bp.
    ))
S4Vectors::mcols(tss)$tss_id <- 1:length(tss)
idx_da_tss <- c(4, 6, 7, 13, 15, 16)
