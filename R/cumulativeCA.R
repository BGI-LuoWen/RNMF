cumulativeCA = function (file = NULL, Filetype = c("MAF", "VCF"), Pfile = NULL, 
    Sfile = NULL, sample.id = "Tumor_Sample_Barcode", chr = "Chromosome", 
    pos = "Start_position", ref = "Reference_Allele", alt = "Tumor_Seq_Allele2", 
    Hugo = "Hugo_Symbol", AnalCOSMICSigType = c("SBS", "DBS", 
        "ID"), genome.build = c("Ch37", "Ch38"), groupFile = NULL, 
    geneListSortFile = NULL, ID.mc.cores = 2, ID.row83 = TRUE, 
    plot = TRUE) 
{
    if (is.null(file)) {
        stop("Please enter a mutation dataset, such as the annotation result file in MAF or VCF format!\n")
    }
    if (is.null(Pfile)) {
        stop("Please enter a mutational signature matrix file P with a format like the output result of RNMF software!\n")
    }
    if (is.null(Sfile)) {
        stop("Please enter a abundance fractions matrix S with a format like the output result of RNMF software!\n")
    }
    library(data.table)
    mutdata = as.data.frame(fread(file))[, c(Hugo, chr, pos, 
        ref, alt, sample.id)]
    ix = grepl("_", mutdata[, chr]) | grepl("M", mutdata[, chr])
    mut = mutdata[!ix, ]
    mut[, chr] = sub("chr", "", mut[, chr])
    mut[, chr] = sub("Chr", "", mut[, chr])
    mut[, chr] = sub("CHR", "", mut[, chr])
    mut[, chr] = paste("chr", mut[, chr], sep = "")
    Nduplicate = mut[, c(Hugo, chr)]
    Nduplicate = unique(Nduplicate[order(Nduplicate[, Hugo]), 
        ])
    Nduplicatenum = t(table(Nduplicate[, Hugo]))
    Nduplicate = colnames(Nduplicatenum)[which(Nduplicatenum > 
        1)]
    index_ix = NULL
    for (i in 1:length(Nduplicate)) {
        ix = which(mut[, Hugo] == Nduplicate[i])
        if (length(ix) > 0) {
            index_ix = c(index_ix, ix)
        }
    }
    if (!is.null(index_ix)) {
        mut = mut[-index_ix, ]
    }
    samMatrixsigsData = samMatrixSigs(Pfile, Sfile, AnalCOSMICSigType)
    geneMatrixNumData = geneMatrixNum(mut, geneListSortFile, 
        Filetype, AnalCOSMICSigType, genome.build, sample.id, 
        chr, pos, ref, alt, Hugo, ID.mc.cores, ID.row83)
    geneList = geneMatrixNumData$geneList
    dir.create("samplesResults", showWarnings = FALSE, recursive = TRUE, 
        mode = "0777")
    Samples = intersect(unlist(samMatrixsigsData$SampleID), unlist(geneMatrixNumData$SampleID))
    n = length(Samples)
    samMatrixsigsDatatmp = list()
    geneMatrixNumDatatmp = list()
    for (i in 1:n) {
        ix = which(unlist(samMatrixsigsData$SampleID) == Samples[i])
        iy = which(unlist(geneMatrixNumData$SampleID) == Samples[i])
        if (length(ix) > 0 && length(iy) > 0) {
            samMatrixsigsDatatmp$SampleID[[i]] = samMatrixsigsData$SampleID[[ix]]
            samMatrixsigsDatatmp$P[[i]] = samMatrixsigsData$P[[ix]]
            geneMatrixNumDatatmp$SampleID[[i]] = geneMatrixNumData$SampleID[[iy]]
            geneMatrixNumDatatmp$samGeneMatrix[[i]] = geneMatrixNumData$samGeneMatrix[[iy]]
        }
    }
    samMatrixsigsData = samMatrixsigsDatatmp
    geneMatrixNumData = geneMatrixNumDatatmp
    rm("geneMatrixNumDatatmp", "samMatrixsigsDatatmp")
    if (is.null(groupFile))
	{
        group = as.data.frame(cbind(Samples, "All"))
    }else{
        library(data.table)
        group = as.data.frame(fread(groupFile))
        index = NULL
        for (i in 1:n) {
            ix = which(group[, 1] == Samples[i])
            if (length(ix) > 0) {
                index = c(index, ix)
            }
        }
        group = group[index, ]
    }
    colnames(group) = c("SampleID", "Group")
    ngs = unique(sort(group$Group))
    ngl = length(ngs)
    asnumeric = function(x) {
        x = as.numeric(as.character(x))
        return(x)
    }
    xy = function(x) {
        a = sum(x)
        if (a == 0) {
            a = 1
        }
        x = x/a
        return(x)
    }
    geneCumulativeContributionAbundance = NULL
    geneRelativeCumulativeContributionAbundance = NULL
    samMatrixDateRes = list()
    id_sam = 1
    geneListtmp = mut[, c(Hugo, chr)]
    data(knownGene)
    geneList1 = unique(merge(geneListtmp, knownGene, by = 1:2, all.x = TRUE, sort = FALSE))
    colnames(geneList1) = c(Hugo, chr, "Length")
    GinNs = intersect(geneList[, Hugo], geneList1[, Hugo])
    allGnames = colnames(geneMatrixNumData$samGeneMatrix[[1]])
    ix = is.na(geneList$Length)
    geneListtmp1 = geneList[!ix, ]
    ix = is.na(geneList1$Length) | grepl("_", geneList1[, chr])
    geneList1tmp1 = geneList1[!ix, ]
    GinNs1 = intersect(geneListtmp1[, Hugo], geneList1tmp1[, Hugo])
    InNs = intersect(geneList1tmp1[, Hugo],allGnames)
    rownames(geneList1tmp1) = geneList1tmp1[, Hugo]
    for (i in 1:ngl)
	{
        dir.create(paste("samplesResults/image_", ngs[i], sep = ""), 
            showWarnings = FALSE, recursive = TRUE, mode = "0777")
        dir.create(paste("samplesResults/image_relative_", ngs[i], 
            sep = ""), showWarnings = FALSE, recursive = TRUE, 
            mode = "0777")
        ix = group$Group == ngs[i]
        ngSamples = group$SampleID[ix]
        for (j in 1:length(ngSamples)) {
            samMatrixDateRes$SampleID[[id_sam]] = ngSamples[j]
            dir.create(paste("samplesResults/", ngSamples[j], 
                "/", AnalCOSMICSigType, sep = ""), showWarnings = FALSE, 
                recursive = TRUE, mode = "0777")
            index = which(samMatrixsigsData$SampleID == ngSamples[j])
            rho = samMatrixsigsData$P[[index]]
            rho = rho[, -1]
            rho = apply(rho, 2, asnumeric)
            index = which(geneMatrixNumData$SampleID == ngSamples[j])
            if (length(index) > 0) {
                G = geneMatrixNumData$samGeneMatrix[[index]]
                G = G[, -1]
                G = apply(G, 2, asnumeric)
                G = t(G)
            }
            else {
                G = matrix(0, nrow(geneList1), nrow(samMatrixsigsData$P[[1]]))
                rownames(G) = geneList1[, Hugo]
            }
            GG = apply(G, 2, xy)
            GG = t(apply(GG, 1, xy))
            theta = GG %*% rho
            theta = theta[GinNs, ]
            tmp = t(apply(theta, 1, xy))
            thetatable = cbind(t(t(rownames(tmp))), tmp)
            colnames(thetatable)[1] = "GeneName"
            write.table(thetatable, file = paste("samplesResults/", 
                ngSamples[j], "/", AnalCOSMICSigType, "/", ngSamples[j], 
                ".geneSigsResult.txt", sep = ""), quote = F, 
                col.names = T, row.names = F, sep = "\t")
            samMatrixDateRes$theta[[id_sam]] = tmp
            if (is.null(geneCumulativeContributionAbundance)) {
                geneCumulativeContributionAbundance = tmp
            }
            else {
                geneCumulativeContributionAbundance = geneCumulativeContributionAbundance + 
                  tmp
            }
            if (plot && !is.null(geneListSortFile)) {
                theta = apply(theta, 1, xy)
                png(filename = paste("samplesResults/image_", 
                  ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                  "_ALL.png", sep = ""), width = 480, height = 480)
                breaks.frequency <- seq(from = min(theta), to = max(theta), 
                  length.out = 10)
                myColors <- colorRampPalette(c("grey98", "red"))
                image(1:nrow(theta), 1:ncol(theta), as.matrix(theta), 
                  breaks = breaks.frequency, col = myColors(length(breaks.frequency) - 
                    1), axes = FALSE, cex = 1.5, xlab = "", ylab = "")
                dev.off()
                nm = length(GinNs)
                kn = nrow(theta)
                Nnameset = rownames(theta)
                for (sigi in 1:kn) {
                  xxda = theta[sigi, ]
                  if (length(xxda) < ceiling(sqrt(nm)) * ceiling(sqrt(nm))) {
                    xxda = c(xxda, rep(0, ceiling(sqrt(nm)) * 
                      ceiling(sqrt(nm)) - length(xxda)))
                  }
                  imageMatrix = matrix(xxda, ceiling(sqrt(nm)), 
                    ceiling(sqrt(nm)), byrow = TRUE)
                  png(filename = paste("samplesResults/image_", 
                    ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                    "_", Nnameset[sigi], ".png", sep = ""), width = 480, 
                    height = 480)
                  breaks.frequency <- seq(from = min(imageMatrix), 
                    to = max(imageMatrix), length.out = 10)
                  myColors <- colorRampPalette(c("grey98", "red"))
                  image(1:nrow(imageMatrix), 1:ncol(imageMatrix), 
                    as.matrix(imageMatrix), breaks = breaks.frequency, 
                    col = myColors(length(breaks.frequency) - 
                      1), axes = FALSE, cex = 1.5, xlab = "", 
                    ylab = "")
                  dev.off()
                }
            }
            G1 = G[InNs, ]
            G1 = G1/geneList1tmp1[InNs, ]$Length
            GG1 = apply(G1, 2, xy)
            GG1 = t(apply(GG1, 1, xy))
            gamma = GG1 %*% rho
            gamma = gamma[GinNs1, ]
            tmp = t(apply(gamma, 1, xy))
            gammatable = cbind(t(t(rownames(tmp))), tmp)
            colnames(gammatable)[1] = "GeneName"
            write.table(gammatable, file = paste("samplesResults/", 
                ngSamples[j], "/", AnalCOSMICSigType, "/", ngSamples[j], 
                ".geneRelativeSigsResult.txt", sep = ""), quote = F, 
                col.names = T, row.names = F, sep = "\t")
            samMatrixDateRes$gamma[[id_sam]] = tmp
            if (is.null(geneRelativeCumulativeContributionAbundance)) {
                geneRelativeCumulativeContributionAbundance = tmp
            }
            else {
                geneRelativeCumulativeContributionAbundance = geneRelativeCumulativeContributionAbundance + 
                  tmp
            }
            if (plot && !is.null(geneListSortFile)) {
                gamma = apply(gamma, 1, xy)
                png(filename = paste("samplesResults/image_relative_", 
                  ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                  "_ALL.png", sep = ""), width = 480, height = 480)
                breaks.frequency <- seq(from = min(gamma), to = max(gamma), 
                  length.out = 10)
                myColors <- colorRampPalette(c("grey98", "red"))
                image(1:nrow(gamma), 1:ncol(gamma), as.matrix(gamma), 
                  breaks = breaks.frequency, col = myColors(length(breaks.frequency) - 
                    1), axes = FALSE, cex = 1.5, xlab = "", ylab = "")
                dev.off()
                nm = length(GinNs1)
                kn = nrow(gamma)
                Nnameset = rownames(gamma)
                for (sigi in 1:kn) {
                  xxda = gamma[sigi, ]
                  if (length(xxda) < ceiling(sqrt(nm)) * ceiling(sqrt(nm))) {
                    xxda = c(xxda, rep(0, ceiling(sqrt(nm)) * 
                      ceiling(sqrt(nm)) - length(xxda)))
                  }
                  imageMatrix = matrix(xxda, ceiling(sqrt(nm)), 
                    ceiling(sqrt(nm)), byrow = TRUE)
                  png(filename = paste("samplesResults/image_relative_", 
                    ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                    "_", Nnameset[sigi], ".png", sep = ""), width = 480, 
                    height = 480)
                  breaks.frequency <- seq(from = min(imageMatrix), 
                    to = max(imageMatrix), length.out = 10)
                  myColors <- colorRampPalette(c("grey98", "red"))
                  image(1:nrow(imageMatrix), 1:ncol(imageMatrix), 
                    as.matrix(imageMatrix), breaks = breaks.frequency, 
                    col = myColors(length(breaks.frequency) - 
                      1), axes = FALSE, cex = 1.5, xlab = "", 
                    ylab = "")
                  dev.off()
                }
            }
            id_sam = id_sam + 1
        }
    }
    geneCumulativeContributionAbundancetable = cbind(t(t(rownames(geneCumulativeContributionAbundance))), 
        geneCumulativeContributionAbundance)
    colnames(geneCumulativeContributionAbundancetable)[1] = "GeneName"
    geneRelativeCumulativeContributionAbundancetable = cbind(t(t(rownames(geneRelativeCumulativeContributionAbundance))), 
        geneRelativeCumulativeContributionAbundance)
    colnames(geneRelativeCumulativeContributionAbundancetable)[1] = "GeneName"
    write.table(geneCumulativeContributionAbundancetable, file = paste("samplesResults/", 
        AnalCOSMICSigType, ".geneCumulativeContributionAbundance.txt", 
        sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
    write.table(geneRelativeCumulativeContributionAbundancetable, 
        file = paste("samplesResults/", AnalCOSMICSigType, ".geneRelativeCumulativeContributionAbundance.txt", 
            sep = ""), quote = F, col.names = T, row.names = F, 
        sep = "\t")
    SampleIDs = unlist(samMatrixDateRes$SampleID)
    nSams = length(SampleIDs)
    nSigs = ncol(samMatrixDateRes$theta[[1]])
    SigIDs = colnames(samMatrixDateRes$theta[[1]])
    geneIDs = rownames(samMatrixDateRes$theta[[1]])
	geneIDs1 = rownames(samMatrixDateRes$gamma[[1]])
    for (i in 1:nSigs)
	{
        dataMatrix0 = NULL
        dataMatrix1 = NULL
        for (j in 1:nSams) {
            dataMatrix0 = cbind(dataMatrix0, t(t(as.numeric(as.character(samMatrixDateRes$theta[[j]][,i])))))
            dataMatrix1 = cbind(dataMatrix1, t(t(as.numeric(as.character(samMatrixDateRes$gamma[[j]][,i])))))
        }
        colnames(dataMatrix0) = SampleIDs
        colnames(dataMatrix1) = SampleIDs
        dataMatrix0 = cbind(t(t(geneIDs)), dataMatrix0)
        colnames(dataMatrix0)[1] = "GeneName"
        write.table(dataMatrix0, file = paste("samplesResults/", 
            SigIDs[i], ".", AnalCOSMICSigType, ".geneCumulativeContributionAbundance.txt", 
            sep = ""), quote = F, col.names = T, row.names = F, 
            sep = "\t")
        dataMatrix1 = cbind(t(t(geneIDs1)), dataMatrix1)
        colnames(dataMatrix1)[1] = "GeneName"
        write.table(dataMatrix1, file = paste("samplesResults/", 
            SigIDs[i], ".", AnalCOSMICSigType, ".geneRelativeCumulativeContributionAbundance.txt", 
            sep = ""), quote = F, col.names = T, row.names = F, 
            sep = "\t")
    }

    if (plot && is.null(geneListSortFile))
	{
        CCAAllFile = paste("samplesResults/", AnalCOSMICSigType, 
            ".geneCumulativeContributionAbundance.txt", sep = "")
        RCCAAllFile = paste("samplesResults/", AnalCOSMICSigType, 
            ".geneRelativeCumulativeContributionAbundance.txt", 
            sep = "")
        CCAAll = as.data.frame(fread(CCAAllFile))
        RCCAAll = as.data.frame(fread(RCCAAllFile))
        CCAAll1 <- CCAAll[, -1]
        CCAAllFileindex <- hclust(dist(CCAAll1), method = "ave")$order
        RCCAAll1 <- RCCAAll[, -1]
        RCCAAllFileindex <- hclust(dist(RCCAAll1), method = "ave")$order
        for (i in 1:nSigs) {
            Stmpfile = paste("samplesResults/", SigIDs[i], ".", 
                AnalCOSMICSigType, ".geneCumulativeContributionAbundance.txt", 
                sep = "")
            RStmpfile = paste("samplesResults/", SigIDs[i], ".", 
                AnalCOSMICSigType, ".geneRelativeCumulativeContributionAbundance.txt", 
                sep = "")
            CCAAll = as.data.frame(fread(Stmpfile))
            RCCAAll = as.data.frame(fread(RStmpfile))
            CCAAll1 <- CCAAll[, -1]
            Stmpfileindex <- hclust(dist(CCAAll1), method = "ave")$order
            RCCAAll1 <- RCCAAll[, -1]
            RStmpfileindex <- hclust(dist(RCCAAll1), method = "ave")$order
            for (i in 1:ngl) {
                ix = group$Group == ngs[i]
                ngSamples = group$SampleID[ix]
                for (j in 1:length(ngSamples)) {
                  SamFile = paste("samplesResults/", ngSamples[j], 
                    "/", AnalCOSMICSigType, "/", ngSamples[j], 
                    ".geneSigsResult.txt", sep = "")
                  CCAAll = as.data.frame(fread(SamFile))
                  theta = CCAAll[CCAAllFileindex, -1]
                  theta = t(theta)
                  png(filename = paste("samplesResults/image_", 
                    ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                    "_ALL.png", sep = ""), width = 480, height = 480)
                  breaks.frequency <- seq(from = min(theta), 
                    to = max(theta), length.out = 10)
                  myColors <- colorRampPalette(c("grey98", "red"))
                  image(1:nrow(theta), 1:ncol(theta), as.matrix(theta), 
                    breaks = breaks.frequency, col = myColors(length(breaks.frequency) - 
                      1), axes = FALSE, cex = 1.5, xlab = "", 
                    ylab = "")
                  dev.off()
                  nm = length(GinNs)
                  kn = nrow(theta)
                  Nnameset = rownames(theta)
                  for (sigi in 1:kn) {
                    xxda = theta[sigi, ]
                    if (length(xxda) < ceiling(sqrt(nm)) * ceiling(sqrt(nm))) {
                      xxda = c(xxda, rep(0, ceiling(sqrt(nm)) * 
                        ceiling(sqrt(nm)) - length(xxda)))
                    }
                    imageMatrix = matrix(xxda, ceiling(sqrt(nm)), 
                      ceiling(sqrt(nm)), byrow = TRUE)
                    png(filename = paste("samplesResults/image_", 
                      ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                      "_", Nnameset[sigi], ".png", sep = ""), 
                      width = 480, height = 480)
                    breaks.frequency <- seq(from = min(imageMatrix), 
                      to = max(imageMatrix), length.out = 10)
                    myColors <- colorRampPalette(c("grey98", 
                      "red"))
                    image(1:nrow(imageMatrix), 1:ncol(imageMatrix), 
                      as.matrix(imageMatrix), breaks = breaks.frequency, 
                      col = myColors(length(breaks.frequency) - 
                        1), axes = FALSE, cex = 1.5, xlab = "", 
                      ylab = "")
                    dev.off()
                  }
                  SamFile = paste("samplesResults/", ngSamples[j], 
                    "/", AnalCOSMICSigType, "/", ngSamples[j], 
                    ".geneRelativeSigsResult.txt", sep = "")
                  CCAAll = as.data.frame(fread(SamFile))
                  gamma = CCAAll[RCCAAllFileindex, -1]
                  gamma = t(gamma)
                  png(filename = paste("samplesResults/image_relative_", 
                    ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                    "_ALL.png", sep = ""), width = 480, height = 480)
                  breaks.frequency <- seq(from = min(gamma), 
                    to = max(gamma), length.out = 10)
                  myColors <- colorRampPalette(c("grey98", "red"))
                  image(1:nrow(gamma), 1:ncol(gamma), as.matrix(gamma), 
                    breaks = breaks.frequency, col = myColors(length(breaks.frequency) - 
                      1), axes = FALSE, cex = 1.5, xlab = "", 
                    ylab = "")
                  dev.off()
                  nm = length(GinNs1)
                  kn = nrow(gamma)
                  Nnameset = rownames(gamma)
                  for (sigi in 1:kn) {
                    xxda = gamma[sigi, ]
                    if (length(xxda) < ceiling(sqrt(nm)) * ceiling(sqrt(nm))) {
                      xxda = c(xxda, rep(0, ceiling(sqrt(nm)) * 
                        ceiling(sqrt(nm)) - length(xxda)))
                    }
                    imageMatrix = matrix(xxda, ceiling(sqrt(nm)), 
                      ceiling(sqrt(nm)), byrow = TRUE)
                    png(filename = paste("samplesResults/image_relative_", 
                      ngs[i], "/", ngSamples[j], "_", AnalCOSMICSigType, 
                      "_", Nnameset[sigi], ".png", sep = ""), 
                      width = 480, height = 480)
                    breaks.frequency <- seq(from = min(imageMatrix), 
                      to = max(imageMatrix), length.out = 10)
                    myColors <- colorRampPalette(c("grey98", 
                      "red"))
                    image(1:nrow(imageMatrix), 1:ncol(imageMatrix), 
                      as.matrix(imageMatrix), breaks = breaks.frequency, 
                      col = myColors(length(breaks.frequency) - 
                        1), axes = FALSE, cex = 1.5, xlab = "", 
                      ylab = "")
                    dev.off()
                  }
                }
            }
        }
    }

}
