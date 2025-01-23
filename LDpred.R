suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(magrittr))
suppressMessages(library(remotes))
suppressMessages(library(bigsnpr))
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
#library(fmsb) #calculate the pseudo R2
option_list = list(
        make_option(c("-f","--fold"), type = "character", default = "/home/wangfan819/PGD/BRI002/gdmPRS/LDpred/", help = "geno and summary statistics direction"),
	make_option(c("-m","--model"), type = "character",default = "infinite",  help = "infinite/grid/auto"),
	make_option(c("-f","--filter"), type = "numeric", default = 0, help = "score top percentage")
)

##get phenos
phenotype = fread("GDM.pheno")
covariate = fread("GDM.cov")
pcs = fread("GDM.0.8.eigenvec")
pcans = ncol(pcs)-2
colnames(pcs) = c("FID","IID", paste0("PC",1:pcans))
pheno = merge(phenotype, covariate) %>% merge(., pcs)

##get summary statistic file
sumstats = bigreadr::fread2("T2D.QC.gz") 
sumstats1 = sumstats
sumstats$BETA = sumstats1$OR
sumstats$OR = sumstats1$BETA
names(sumstats) = c("chr","pos","rsid","a1","a0","n_eff","beta_se","p","OR","INFO","MAF","beta")
nonFlipN = dim(sumstats)[1]/2
sumstats = sumstats[1:nonFlipN,]
sumstats = sumstats[sumstats$chr != "X",]
sumstats$chr = as.integer(sumstats$chr)
#str(sumstats)

if(FALSE){
##get HapMap3 SNPs
info = readRDS("/home/wangfan819/PGD/BRI002/gdmPRS/LDpred/map_hm3_ldpred2.rds")
#info = readRDS(runonce::download_file("https://ndownloader.figshare.com/files/25503788",fname = "map_hm3_ldpred2.rds"))
info$SNPid = paste(info$chr,info$pos,info$a1,info$a0,sep=":")
info$SNPid2 = paste(info$chr,info$pos,info$a0,info$a1,sep=":")#Flip A1 A2 to cover the chip SNP
sumstats = sumstats[(sumstats$rsid %in% info$SNPid) | (sumstats$rsid %in% info$SNPid2),]
}

##Calculate the LD matrix
NCORES = nb_cores()#Get maximum amount of cores
tmp = tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)#Open a temporary file
corr = NULL
ld = NULL#Initialize variables for storing the LD score and LD matrix
fam.order = NULL
snp_readBed("GDM.0.8.QC.bed")#preprocess the bed file (only need to do once for each data set)
obj.bigSNP = snp_attach("GDM.0.8.QC.rds")#attach the genotype object
#str(obj.bigSNP, max.level = 2, strict.width = "cut")
map = obj.bigSNP$map[-3]
names(map) = c("chr", "rsid", "pos", "a1", "a0")#extract the SNP information from the genotype
#sumstats$chr = as.numeric(sumstats$chr)#no X

info_snp = snp_match(sumstats, map, join_by_pos = FALSE,strand_flip =FALSE, match.min.prop = 0.1)#map base and target, remove_dups default TRUE
#info_snp1 = snp_match(sumstats, map,strand_flip =FALSE,remove_dups = FALSE,match.min.prop = 0.1)

genotype1 = obj.bigSNP$genotypes#Assign the genotype to a variable for easier downstream analysis
genotype = snp_fastImputeSimple(genotype1)
CHR = map$chr
POS = map$pos#Rename the data structures
POS2 = snp_asGeneticPos(CHR, POS, dir = "./1000genomesCM/")#download the 1000G file to the current directory (".")
#POS2 = obj.bigSNP$map$genetic.dist


# calculate LD
tmp = tempfile(tmpdir = "tmp-data")
for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr = which(info_snp$chr == chr)
    ind.chr2 = info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD,length = chrAsnps^2
    corr0 = snp_cor(genotype, ind.col = ind.chr2, ncores = NCORES, infos.pos = POS2[ind.chr2], size = 3 / 1000)
    if (chr == 1) {
        ld = Matrix::colSums(corr0^2)
        corr = as_SFBM(corr0, tmp)#corr = as_SFBM(corr0, tmp, compact = TRUE)
    } else {
        ld = c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))#corr: N*N, N=number of snp in info_snp
    }
}
fam.order = as.data.table(obj.bigSNP$fam)
setnames(fam.order,c("family.ID", "sample.ID"),c("FID", "IID"))
##LD score regression
df_beta = info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc = snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est = ldsc[["h2"]]#ldsc with two columns: int  h2

##Calculate the null R2
y = pheno[fam.order, on = c("FID", "IID")]#y is of the same order as the sample ordering in the genotype file
null.model = paste("PC", 1:pcans, sep = "", collapse = "+") %>%
    paste0("T2D~Sex+", .) %>%
    as.formula %>%
    glm(., data = y, family=binomial) %>%
    summary
unVar1 = ncol(y); unVar2 = ncol(y) -3
y3 = as.data.frame(y[,-(unVar2:unVar1)])
null.model1 = glm(T2D~1,family = binomial(),data=y3[,!colnames(y3)%in%c("FID","IID")])
null.model2 = glm(T2D~.,family = binomial(),data=y3[,!colnames(y3)%in%c("FID","IID")])
#null.r2 = fmsb::NagelkerkeR2(null.model)  #Nagelkerke R2 is biased when there are ascertainment of samples
null.r2 = as.numeric(1-logLik(null.model2)/logLik(null.model1))#McFadden’s Pseudo-R2

###Obtain model PRS
##infinitesimal model
beta_inf = snp_ldpred2_inf(corr, df_beta, h2 = h2_est)#numeric, length = map SNPs(bed file)
if(is.null(obj.bigSNP)){
    obj.bigSNP = snp_attach("GDM.0.8.QC.rds")
}
genotype = obj.bigSNP$genotypes #row:samples ; col:snps of map
ind.test = 1:nrow(genotype) #sample number
pred_inf = big_prodVec( genotype , beta_inf, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)

##grid model
p_seq = signif(seq_log(1e-4, 1, length.out = 17), 2)
h2_seq = round(h2_est * c(0.7, 1, 1.4), 4)
grid.param = expand.grid(p = p_seq,h2 = h2_seq,sparse = c(FALSE, TRUE))
beta_grid = snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
if(is.null(obj.bigSNP)){
    obj.bigSNP = snp_attach("GDM.0.8.QC.rds")
}
genotype = obj.bigSNP$genotypes
ind.test = 1:nrow(genotype)
pred_grid = big_prodMat( genotype , beta_grid, ind.col = info_snp$`_NUM_ID_`)

##auto model


#####Get the final performance of the LDpred models
reg.formula = paste("PC", 1:pcans, sep = "", collapse = "+") %>%
    paste0("T2D~PRS+Sex+", .) %>%
    as.formula
reg.dat = y
reg.dat$PRS = pred_inf
inf.model = glm(reg.formula, dat=reg.dat) %>%
    summary


SNPmerge = function (sumstats, info_snp, strand_flip = TRUE, join_by_pos = TRUE,
    remove_dups = TRUE, match.min.prop = 0.5)
{
    sumstats <- as.data.frame(sumstats)
    info_snp <- as.data.frame(info_snp)
    sumstats$`_NUM_ID_` <- rows_along(sumstats)
    info_snp$`_NUM_ID_` <- rows_along(info_snp)
    min_match <- match.min.prop * min(nrow(sumstats), nrow(info_snp))
    join_by <- c("chr", NA, "a0", "a1")
    join_by[2] <- if (join_by_pos)
        "pos"
    else "rsid"
    if (!all(c(join_by, "beta") %in% names(sumstats)))
        stop2("Please use proper names for variables in 'sumstats'. Expected '%s'.",
            paste(c(join_by, "beta"), collapse = ", "))
    if (!all(c(join_by, "pos") %in% names(info_snp)))
        stop2("Please use proper names for variables in 'info_snp'. Expected '%s'.",
            paste(unique(c(join_by, "pos")), collapse = ", "))
    message2("%s variants to be matched.", format(nrow(sumstats),
        big.mark = ","))
    sumstats <- sumstats[vctrs::vec_in(sumstats[, join_by[1:2]],
        info_snp[, join_by[1:2]]), ]
    if (nrow(sumstats) == 0)
        stop2("No variant has been matched.")
    if (strand_flip) {
        is_ambiguous <- with(sumstats, paste(a0, a1) %in% c("A T",
            "T A", "C G", "G C"))
        message2("%s ambiguous SNPs have been removed.", format(sum(is_ambiguous),
            big.mark = ","))
        sumstats2 <- sumstats[!is_ambiguous, ]
        sumstats3 <- sumstats2
        sumstats2$`_FLIP_` <- FALSE
        sumstats3$`_FLIP_` <- TRUE
        sumstats3$a0 <- flip_strand(sumstats2$a0)
        sumstats3$a1 <- flip_strand(sumstats2$a1)
        sumstats3 <- rbind(sumstats2, sumstats3)
    }
    else {
        sumstats3 <- sumstats
        sumstats3$`_FLIP_` <- FALSE
    }
    sumstats4 <- sumstats3
    sumstats3$`_REV_` <- FALSE
    sumstats4$`_REV_` <- TRUE
    sumstats4$a0 <- sumstats3$a1
    sumstats4$a1 <- sumstats3$a0
    sumstats4$beta <- -sumstats3$beta
    sumstats4 <- rbind(sumstats3, sumstats4)
    matched <- merge(as.data.table(sumstats4), as.data.table(info_snp),
        by = join_by, all = FALSE, suffixes = c(".ss", ""))
    if (remove_dups) {
        dups <- vctrs::vec_duplicate_detect(matched[, c("chr",
            "pos")])
        if (any(dups)) {
            matched <- matched[!dups, ]
            message2("Some duplicates were removed.")
        }
    }
    message2("%s variants have been matched; %s were flipped and %s were reversed.",
        format(nrow(matched), big.mark = ","), format(sum(matched$`_FLIP_`),
            big.mark = ","), format(sum(matched$`_REV_`), big.mark = ","))
    if (nrow(matched) < min_match)
        stop2("Not enough variants have been matched.")
    as.data.frame(matched[, `:=`(c("_FLIP_", "_REV_"), NULL)][order(chr,
        pos)])
}

CMinfo = function (infos.chr, infos.pos, dir = tempdir(), ncores = 1,
    rsid = NULL)
{
    assert_package("R.utils")
    assert_lengths(infos.chr, infos.pos)
    if (!is.null(rsid))
        assert_lengths(rsid, infos.pos)
    snp_split(infos.chr, function(ind.chr, pos, dir, rsid) {
        chr <- attr(ind.chr, "chr")
        basename <- paste0("chr", chr, ".OMNI.interpolated_genetic_map")
        mapfile <- file.path(dir, basename)
        if (!file.exists(mapfile)) {
            url <- paste0("https://github.com/joepickrell/1000-genomes-genetic-maps/",
                "raw/master/interpolated_OMNI/", basename, ".gz")
            gzfile <- paste0(mapfile, ".gz")
            utils::download.file(url, destfile = gzfile, quiet = TRUE)
            R.utils::gunzip(gzfile)
        }
        map.chr <- bigreadr::fread2(mapfile, showProgress = FALSE)
        if (is.null(rsid)) {
            ind <- bigutilsr::knn_parallel(as.matrix(map.chr$V2),
                as.matrix(pos[ind.chr]), k = 1, ncores = 1)$nn.idx
            new_pos <- map.chr$V3[ind]
        }
        else {
            ind <- match(rsid[ind.chr], map.chr$V1)
            new_pos <- map.chr$V3[ind]
            indNA <- which(is.na(ind))
            if (length(indNA) > 0) {
                pos.chr <- pos[ind.chr]
                new_pos[indNA] <- suppressWarnings(stats::spline(pos.chr,
                  new_pos, xout = pos.chr[indNA], method = "hyman")$y)
            }
        }
        new_pos
    }, combine = "c", pos = infos.pos, dir = dir, rsid = rsid,
        ncores = ncores)
}

big_prodVec = function (X, y.col, ind.row = rows_along(X), ind.col = cols_along(X),
    center = NULL, scale = NULL, ncores = 1)
{
    check_args()
    assert_lengths(y.col, ind.col)
    if (length(ind.row) == 0 || length(ind.col) == 0)
        return(rep(0, length(ind.row)))
    if (!is.null(scale)) {
        assert_lengths(scale, ind.col)
        y.col <- y.col/as_vec(scale)
    }
    if (!is.null(center)) {
        assert_lengths(center, ind.col)
        center2 <- drop(crossprod(as_vec(center), y.col))
    }
    res <- pMatVec4(X, y.col, ind.row, ind.col, ncores = ncores)
    if (is.null(center))
        res
    else res - center2
}
