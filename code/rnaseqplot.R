#Initial step requires STAR-HTseq running
#read in files using tximport, replace path with local path
library(DESeq2)
library(tximportData)
library(GenomicFeatures)
library(ggplot2)

ff <-
  list.files(path = "./CasMINI/RNAseq_dscGFP/",
             pattern = "*ReadsPerGene.out.tab$",
             full.names = TRUE)
counts.files <- lapply(ff, read.table, skip = 4)

#store counts from tab files
#use column two for strand independent data, 3 for 1st strand, 4 for reverse strand
counts <-
  as.data.frame(sapply(counts.files, function(x)
    x[, 2]))
ff <- gsub("ReadsPerGene[.]out[.]tab", "", ff)
ff <- gsub("[.]/CasMINI/RNAseq_dscGFP/", "", ff)
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1

#call DESeq2 and invoke
coldataCasMINI <-
  data.frame(
    group = c(rep("dCasMINI", 4), rep("dLbCas12", 4)),
    treatment = rep(c(rep("sglacZ", 2), rep("sgTet", 2)), 2),
    replicate = rep(rep(c(1, 2), 2), 2),
    mode = rep("PE", 8),
    row.names = colnames(counts)
  )
coldataCasMINI$all <-
  as.factor(paste(coldataCasMINI$group, coldataCasMINI$treatment, sep = "_"))
coldataCasMINI$x <- as.factor(c("B", "B", "C", "C", "F", "F", "G", "G"))

ddsCasMINI <-
  DESeqDataSetFromMatrix(counts, colData = coldataCasMINI, design = ~ x)
ddsCasMINI <- DESeq(ddsCasMINI)
countsCasMINI <- as.data.frame(counts(ddsCasMINI, normalize = TRUE))

#function for linear regression (statistics) and coefficients of determination
lm_eqn <- function(df) {
  m <- lm(y ~ x, df)
  
  eq <-
    substitute(
      italic(y) == a + b %.% italic(x) * "," ~  ~ italic(r) ^ 2 ~ "=" ~ r2,
      list(
        a = format(unname(coef(m)[1]), digits = 2),
        b = format(unname(coef(m)[2]), digits = 2),
        r2 = format(summary(m)$r.squared, digits = 3)
      )
    )
  as.character(as.expression(eq))
  
}
Rsquaredval <- function(df) {
  colnames(df) <- c("x", "y")
  summary(lm(y ~ x, df))$r.squared
}

#Pearson's correlatiom
cor.test(data.frame(countsCasMINI[, 1:2]), method = c("pearson"))

#plot scatter graphs (with statistics)
scatter_ab <-
  ggplot(log10(countsCasMINI), aes(x = XX26_1, y = XX26_2)) + geom_point(shape =
                                                                           1) + geom_smooth(method = lm)
scatter_ab2 <-
  scatter_ab + geom_text(
    x = 25,
    y = 300,
    label = lm_eqn(countsCasMINI[, 1:2]),
    parse = TRUE
  )

#plotscatter function for making a scatter plot with linear regression
plotscatter <- function(df) {
  colnames(df) <- c("x", "y")
  p <- ggplot(data = log10(df + 1), aes(x = x, y = y)) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black",
      formula = y ~ x
    ) +
    geom_point()
  p1 <- p + geom_text(
    x = 1,
    y = 3,
    label = lm_eqn(df),
    parse = TRUE
  )
  p1
}
#plotscatter run function
plotscatter(countsCasMINI[, c(1, 2)])

rsquaredvals <-
  cbind(sapply(1:8, function(x)
    sapply(1:8, function(y)
      Rsquaredval(countsCasMINI[, c(x, y)]))))
colnames(rsquaredvals) <-
  rownames(rsquaredvals) <- colnames(countsCasMINI)

#Colored plotting of specific genes example
plot(
  log10(countsCasMINI[, c(1, 3)]),
  col = ifelse(countsCasMINI2$gene_name == "EGFP", "green", "black"),
  pch = 20
)

#tpm conversion function
tpm3 <- function(counts, len) {
  x <- counts / len
  return(t(t(x) * 1e6 / colSums(x)))
}
#This code generates tpm normalization.
txdb <-
  makeTxDbFromGFF("gencode.v37.chr_patch_hapl_scaff.annotation.withdscGFP.gtf",
                  format = "auto")

exons.list.per.gene <- exonsBy(txdb, by = "gene")
write.table(exons.list.per.gene, "exons.list.per.gene")
exons.list.per.gene <- read.table("exons.list.per.gene", header = TRUE)

exonic.gene.sizes <-
  as.data.frame(sum(width(reduce(
    exons.list.per.gene
  ))))

CasMINItpm <-
  tpm3(countsCasMINI, exonic.gene.sizes[rownames(countsCasMINI), 1])

CasMINItpmmerged <-
  10 ^ cbind(
    apply(log10(CasMINItpm + 0.01)[, 1:2], 1, mean),
    apply(log10(CasMINItpm + 0.01)[, 3:4], 1, mean),
    apply(log10(CasMINItpm + 0.01)[, 5:6], 1, mean),
    apply(log10(CasMINItpm + 0.01)[, 7:8], 1, mean)
  ) - 0.01

#simple scatter plots
plot(
  log10(CasMINItpm[, c(1, 3)] + 0.01),
  col = ifelse(countsCasMINI2$gene_name == "EGFP", "red3", "red3"),
  pch = 20,
  cex = ifelse(countsCasMINI2$gene_name == intgene, 1, 0.25),
  xlim = c(1, 4),
  ylim = c(1, 4)
)
points(
  log10(CasMINItpm[, c(5, 7)] + 0.01),
  col = ifelse(countsCasMINI2$gene_name == "EGFP", "slateblue1", "slateblue1"),
  pch = 20,
  cex = ifelse(countsCasMINI2$gene_name == intgene, 1, 0.25),
  xlim = c(1, 4),
  ylim = c(1, 4)
)

#Optional: obtain real gene names from gencode and change thenames
annotations <-
  read.table(
    "./CasMINI/RNAseq_dscGFP/gencode.v37.chr_patch_hapl_scaff.annotation.withdscGFP.gtf",
    sep = "\t"
  )

library(stringr)
genenames2 <-
  stringr::str_split_fixed(annotations[, 9], pattern = ";", n = 15)
annogenenames <-
  apply(genenames2, 1, function(x)
    str_split_fixed(x[grepl("gene_name", x)], pattern = " ", n = 3)[3])
annogeneids <-
  apply(genenames2, 1, function(x)
    str_split_fixed(x[grepl("gene_id", x)], pattern = " ", n = 3)[2])
annomasterkey <- unique(cbind(annogeneids, annogenenames))

rownames(annomasterkey) <- annomasterkey[, 1]
countsCasMINI2 <- countsCasMINI
countsCasMINI2$gene_name <-
  annomasterkey[rownames(countsCasMINI2), 2]

#This part generates Scatterplots with gene label.
intgene <- "EGFP"
intplot <-
  as.data.frame(rbind(log10(CasMINItpmmerged[, c(1, 2)] + 1), log10(CasMINItpmmerged[, c(3, 4)] +
                                                                      1)))
colnames(intplot) <- c("sgNT", "sgT")
intplot$color <-
  c(rep("red3", nrow(CasMINItpmmerged)), rep("gray", nrow(CasMINItpmmerged)))
intplot$cex <- ifelse(countsCasMINI2$gene_name == "EGFP", 2, 0.5)
intplot2 <- intplot[sample(nrow(intplot)), ]
plot(
  intplot2[, c("sgNT", "sgT")],
  col = alpha(intplot2$color, 0.5),
  pch = 20,
  cex = intplot2$cex,
  xlim = c(1, 5),
  ylim = c(1, 5),
  xlab = c("sgLac"),
  ylab = c("sgGFP")
)

#final plot function with colored gene
plot(
  intplot2[, c("Mini", "Cas12")],
  col = alpha(intplot2$color, 0.5),
  pch = 20,
  cex = intplot2$cex,
  xlim = c(1, 5),
  ylim = c(1, 5),
  xlab = c("Mini"),
  ylab = c("Cas12")
)


#violin plot functions
GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    # Original function by Jan Gleixner (@jan-glx)
    # Adjustments by Wouter van der Bijl (@Axeman)
    data <-
      transform(
        data,
        xminv = x - violinwidth * (x - xmin),
        xmaxv = x + violinwidth * (xmax - x)
      )
    grp <- data[1, "group"]
    newdata <-
      plyr::arrange(transform(data, x = if (grp %% 2 == 1)
        xminv
        else
          xmaxv), if (grp %% 2 == 1)
            y
        else-y)
    newdata <-
      rbind(newdata[1,], newdata, newdata[nrow(newdata),], newdata[1,])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
      round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 &
        !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <-
        create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
      aesthetics <-
        data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <-
        rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <-
        GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

create_quantile_segment_frame <-
  function(data,
           draw_quantiles,
           split = FALSE,
           grp = NULL) {
    dens <- cumsum(data$density) / sum(data$density)
    ecdf <- stats::approxfun(dens, data$y)
    ys <- ecdf(draw_quantiles)
    violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
    violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
    violin.xs <- (stats::approxfun(data$y, data$x))(ys)
    if (grp %% 2 == 0) {
      data.frame(
        x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
        y = rep(ys, each = 2),
        group = rep(ys, each = 2)
      )
    } else {
      data.frame(
        x = ggplot2:::interleave(violin.xminvs, violin.xs),
        y = rep(ys, each = 2),
        group = rep(ys, each = 2)
      )
    }
  }

geom_split_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "identity",
           ...,
           draw_quantiles = NULL,
           trim = TRUE,
           scale = "area",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomSplitViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        draw_quantiles = draw_quantiles,
        na.rm = na.rm,
        ...
      )
    )
  }

#collect data and draw violin plots
my_data1 <-
  c(apply(CasMINItpm, 1, function(x)
    ifelse(sum(x[c(1, 2, 3, 4)]) <= 5, -1, sd(log2(
      x[c(1, 2, 3, 4)] + 1
    )))))
my_data2 <-
  c(apply(CasMINItpm, 1, function(x)
    ifelse(sum(x[c(5, 6, 7, 8)]) <= 5, -1, sd(log2(
      x[c(5, 6, 7, 8)] + 1
    )))))
my_dataxx <-
  data.frame(as.numeric(cbind(c(my_data1, my_data2))), c(rep("Mini", length(my_data1)), rep("Cas12", length(my_data2))))

colnames(my_dataxx) <- c("y", "x")
my_dataxx$m <- my_dataxx$x
my_dataxx$x <- rep("gfpvsnt", nrow(my_dataxx))
head(my_data)

my_dataxx_mini <- my_dataxx[my_dataxx$m == "Mini" &
                              my_dataxx$y >= 0, ]
my_dataxx_cas12 <-
  my_dataxx[my_dataxx$m == "Cas12" & my_dataxx$y >= 0, ]
my_dataxx4 <-
  rbind(my_dataxx_mini[my_dataxx_mini[, 1] > quantile(my_dataxx_mini[, 1], 0.05) &
                         my_dataxx_mini[, 1] < quantile(my_dataxx_mini[, 1], 0.95), ],
        my_dataxx_cas12[my_dataxx_cas12[, 1] > quantile(my_dataxx_cas12[, 1], 0.05) &
                          my_dataxx_cas12[, 1] < quantile(my_dataxx_cas12[, 1], 0.95), ])

my_dataxx4 <- rbind(my_dataxx_mini[my_dataxx_mini[, 1] < 3, ],
                    my_dataxx_cas12[my_dataxx_cas12[, 1] < 3, ])

ggplot(my_dataxx4[sample(1:nrow(my_dataxx4), 100), ], aes(x, y, fill = m)) + geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))
