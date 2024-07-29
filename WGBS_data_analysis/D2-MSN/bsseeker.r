#!/usr/bin/Rscript


library(MethylSeekR)
        set.seed(123)
        library("BSgenome.Mmusculus.UCSC.mm10")

    for (spl in c("EM_D2_merge.CpG.ss.cov", "WGBS_D2_merge.CpG.ss.cov" )){

        sLengths = seqlengths(Mmusculus)
        sLengths
        meth.gr = readMethylome(FileName=paste0(spl,".bsseekr"),seqLengths=sLengths)

        snps.gr = readSNPTable(FileName="/home/yli/ref/bed/mouse.common.snp", seqLengths=sLengths)
        meth.gr <- removeSNPs(meth.gr, snps.gr)

#PMD segmentation
        svg(paste0(spl,"chr1_distribution.svg"))
        plotAlphaDistributionOneChr(m=meth.gr, chr.sel="chr1", num.cores=1)
        dev.off()
        PMDsegments.gr <- segmentPMDs(m=meth.gr, chr.sel="chr19", seqLengths=sLengths, num.cores=1)

        svg(paste0(spl,"chr19_dis.svg"))
        plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.gr, PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), chr.sel="chr19", num.cores=1)
        dev.off()
        svg(paste0(spl,"PMDSeg.svg"))
        plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr)
        dev.off()

        savePMDSegments(PMDs=PMDsegments.gr,GRangesFilename=paste0(spl,"PMDs.gr.rds"), TableFilename=paste0(spl,"PMDs.tab"))

#UMR and LMR

        library(rtracklayer)
        session <- browserSession()
        genome(session) <- "mm10"
        query <- ucscTableQuery(session, "cpgIslandExt")
        CpGislands.gr <- track(query)
        genome(CpGislands.gr) <- NA
        CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))
        stats <- calculateFDRs(m=meth.gr, CGIs=CpGislands.gr, PMDs=PMDsegments.gr, num.cores=1)

        FDR.cutoff <- 5
        m.sel <- 0.5
        n.sel=as.integer(names(stats$FDRs[as.character(m.sel), ][stats$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
        n.sel

        UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, PMDs=PMDsegments.gr, num.cores=1, myGenomeSeq=Mmusculus, seqLengths=sLengths)
        head(UMRLMRsegments.gr)
        svg(paste0(spl,"final_seg.svg"))
        plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr, PMDs=PMDsegments.gr,meth.cutoff=m.sel)
        dev.off()
        saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename=paste0(spl,"UMRsLMRs.gr.rds"), TableFilename=paste0(spl,"UMRsLMRs.tab"))
    }
